! $Id:$
!
!   This module reads ahi files
!
!

module cx_read_ahi_mod
   
   type, public :: ahi_config_type
      character ( len = 255 ) :: file_base
      character ( len = 400 ) :: data_path
      logical :: chan_on(16)
      character (len = 500) :: filename ( 16)
      character ( len =20) :: varname (16)
      integer :: h5_offset(2)
      integer :: h5_count(2)
   end type ahi_config_type
   
   type ahi_chn_type
      logical :: is_read
      logical :: is_solar_channel
      real,dimension (:,:) , allocatable :: ref
      real,dimension (:,:) , allocatable :: rad
      real,dimension (:,:) , allocatable :: bt
   end type ahi_chn_type
   
   type :: geo_str
      real , dimension (:,:) , allocatable :: solzen
      real , dimension (:,:) , allocatable :: satzen
      real , dimension (:,:) , allocatable :: solaz
      real , dimension (:,:) , allocatable :: sataz  
      real , dimension (:,:) , allocatable :: relaz 
      real , dimension (:,:) , allocatable :: lat
      real , dimension (:,:) , allocatable :: lon  
      real , dimension (:,:) , allocatable :: glintzen
      real , dimension (:,:) , allocatable :: scatangle
      real , dimension (:)   , allocatable :: scan_time
      logical, dimension (:,:), allocatable :: is_space
   end type  geo_str
   
   type, public :: ahi_data_out_type
      type ( ahi_chn_type ) :: chn(16)
      type ( geo_str ) :: geo
      
      contains
      procedure :: deallocate_all
   end type ahi_data_out_type
   
  
  !             
contains

   subroutine get_ahi_data ( config , out )
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: out
      
      call check_input ( config )
      call read_ahi_level1b ( config , out )
   
   end subroutine get_ahi_data
   !
   !
   !
   subroutine check_input ( config)
      use string_functions, only: replace_text
      
      type ( ahi_config_type ) :: config
      
      integer :: i
      character(len=2) :: identifier
      character ( len=255) :: file_for_this_channel

      do i = 1 , 16
        
         write (identifier , fmt ='(i2.2)') i
         
         file_for_this_channel = replace_text ( config % file_base,'B01','B'//identifier)
         config % filename ( i ) = trim (config % data_path)//trim(file_for_this_channel) 
         config % varname ( i) = '/RAD'
         
      end do
      
  
   
   end subroutine 
   
   ! ------------------------------------------------
   !
   ! -----------------------------------------------
   subroutine read_ahi_level1b ( config, ahi )
      !use geo_sat_navigation_mod
      !use hdf5
      use readh5dataset
      use viewing_geometry_module
      
      implicit none
      
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: ahi
      real :: lat ,lon
      integer :: status
      real , pointer :: rad(:,:)          
      integer(kind = 2), dimension ( : , : ) , pointer :: i2d_buffer
      integer:: i_chn
      character (len=120) :: sds_name , attr_name
      real (8) :: scale_factor, add_offset
      integer ( kind = 2 ) :: fillvalue
      integer ( kind = 4 ) , allocatable :: buffer_fake_i4 (:,:)
     
      real ( 8 ) :: cprime 
      logical :: is_solar_channel = .false.
         real(8) :: xx , yy
   real (8) :: cfac
   real (8):: coff
   real (8) :: lfac
   real (8) :: loff
   real(8) :: sub_lon , sub_lat  
   real(8) :: latx, lonx
   real(8) :: fargc
   integer :: ii,jj
   
   real :: GEO_ALTITUDE = 35786.0 !km
      ! - executable
      
      ! - navigation
      
       xx=2300.
   yy=2500.
   cfac = 20466276.
   coff = 2750.
   lfac = 20466276.
   loff = 2750.
   sub_lon =  140.70
   sub_lat = 0.
      
      allocate (  ahi % geo % lon  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % lat  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % solzen  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % solaz  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % satzen  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % sataz  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % relaz  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % glintzen  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % scatangle  (config % h5_count(1),config % h5_count(2)) )
      allocate (  ahi % geo % is_space (config % h5_count(1),config % h5_count(2)) )
      ahi % geo % lon          =   -999.
      ahi % geo % lat      =   -999.
      ahi % geo % solzen      =   -999.
      ahi % geo % solaz      =   -999.
      ahi % geo % satzen      =   -999.
      ahi % geo % sataz      =   -999.
      ahi % geo % relaz      =   -999.
      ahi % geo % glintzen     =   -999.
      ahi % geo % scatangle    =   -999.
      ahi % geo % is_space = .true.
      
      
      print*,config %h5_count,config %h5_offset
      
      do jj =+1 , config % h5_count(2)
     
         do ii =1 ,   config % h5_count(1)
           
            call fgf_to_earth ( 3, dble(ii), dble(jj) , cfac, coff, lfac, loff, sub_lon &
               , lonx , latx )
           
            ahi % geo % lat (ii,jj) = latx  
            ahi % geo % lon (ii,jj) = lonx  
            
            if ( lonx == -999. ) cycle
            
            call  possol (123 ,  2.4  , real(lonx) , real(latx),ahi % geo % solzen (ii,jj),ahi % geo % solaz (ii,jj) )
            
            
            ahi % geo % satzen (ii,jj) = sensor_zenith ( GEO_ALTITUDE, real(sub_lon),real(sub_lat),real(lonx) , real(latx) )
             
            ahi % geo % sataz (ii,jj) = sensor_azimuth (  real(sub_lon),real(sub_lat),real(lonx) , real(latx) )
            
            ahi % geo % relaz(ii,jj) = relative_azimuth (ahi % geo % solaz (ii,jj) , ahi % geo % sataz (ii,jj))
           
            ahi % geo % glintzen(ii,jj) = glint_angle (ahi % geo % solzen (ii,jj) ,ahi % geo % satzen (ii,jj) &
                  ,  ahi % geo % relaz(ii,jj) )
           
            ahi % geo % scatangle(ii,jj) = scattering_angle ( ahi % geo % solzen (ii,jj)  &
                  , ahi % geo % satzen (ii,jj), ahi % geo % relaz(ii,jj))
  
         end do
      end do   
      
      where ( ahi % geo % lat .gt. -299.0)
         ahi % geo % is_space = .false.
      end where
         
   
     
      
      ! - channel data read 
      do i_chn = 1 ,16
       
         if ( .not. config % chan_on ( i_chn ) ) cycle
         
         
         !print*,'Read in AHI FIle > ', trim(config % filename ( i_chn ))
         call h5readdataset ( trim(config % filename ( i_chn ) ) , trim ( config % varname(i_chn) ) &
               , config % h5_offset,config % h5_count, i2d_buffer )
         allocate ( buffer_fake_i4 (config % h5_count(1),config % h5_count(2)))
         
               ! - fortran does not support unsigned integer
         buffer_fake_i4 = i2d_buffer
         where ( buffer_fake_i4 < 0 )
            buffer_fake_i4 = buffer_fake_i4 + 65536
         end where
         
         attr_name = trim(config % varname(i_chn))//'/scale_factor'
         call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), scale_factor )
         attr_name = trim(config % varname(i_chn))//'/add_offset'
         call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), add_offset )
         attr_name = trim(config % varname(i_chn))//'/_FillValue'
         call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), fillvalue )
         if ( fillvalue < 0 ) fillvalue = fillvalue + 65536
         allocate ( ahi % chn(i_chn) % rad (config % h5_count(1),config % h5_count(2)))
          
          ahi % chn(i_chn) % rad = (buffer_fake_i4 * scale_factor) + add_offset
         where ( buffer_fake_i4 == fillvalue )
            ahi % chn(i_chn) % rad = -999.
         end where
         if (allocated ( buffer_fake_i4 ) )  deallocate ( buffer_fake_i4 )
         
         is_solar_channel = .false.
         if ( i_chn < 7 ) is_solar_channel = .true.
         if ( is_solar_channel ) then
            attr_name = trim(config % varname(i_chn))//'/cprime'
            call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), cprime )
            allocate ( ahi % chn(i_chn) % ref (config % h5_count(1),config % h5_count(2)))
            ahi % chn(i_chn) % ref =100. *  ahi % chn(i_chn) % rad * cprime
         end if
        
         ahi % chn(i_chn) % is_read = .true.
         i2d_buffer => null()
         
      end do
 
   end subroutine read_ahi_level1b
   
   
   subroutine deallocate_all (this )
      class ( ahi_data_out_type ) :: this
      integer :: i_chn
      if (allocated ( this % geo % lon)) deallocate ( this % geo % lon) 
      if (allocated ( this % geo % lat)) deallocate ( this % geo % lat) 
      
           
      deallocate (  this % geo % solzen   )
      deallocate (  this % geo % solaz   )
      deallocate (  this % geo % satzen  )
      deallocate (  this % geo % sataz  )
      deallocate (  this % geo % relaz   )
      deallocate (  this % geo % glintzen   )
      deallocate (  this % geo % scatangle   )
      deallocate (  this % geo % is_space  )
      
      do i_chn = 1 ,16
         if (allocated (  this  % chn(i_chn) % rad )) deallocate (this  % chn(i_chn) % rad)
         if (allocated (  this  % chn(i_chn) % bt ) ) deallocate (this  % chn(i_chn) % bt)
         if (allocated (  this  % chn(i_chn) % ref )) deallocate (this  % chn(i_chn) % ref)
      end do
   
   
   end subroutine deallocate_all

end module cx_read_ahi_mod
