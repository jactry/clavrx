! $Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: cx_read_ahi_mod.f90 ( module)
!       viirs_clavrx_bridge (program)
!
! PURPOSE: AHI Reader for NCDF4 / HDF5 files generated at CIMSS 
!
! DESCRIPTION: 
!
!  DEPENDENCIES: 
!     MODULES:
!           readhdf5dataset
!           string_functions
!           viewing_geometry_module
!     SUBROUTINES
!      fgf_to_earth ( c-routine )
!
!
!
! AUTHORS:
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
! REVISION HISTORY:   Created 6 Feb 2015 (AW )
!
! NOTES:
!  AHI  MODIS(CLAVRX)  Wavelength
!           mapping
!    M1    -     8    -    0.412
!    M2    -     9    -    0.445
!    M3    -     3    -    0.488
!    M4    -     4    -    0.555
!    M5    -     1    -    0.672
!    M6    -    15    -    0.746
!    M7    -     2    -    0.865
!    M8    -     5    -    1.240
!    M9    -    26    -    1.378
!    M10   -     6    -    1.610
!    M11   -     7    -    2.250
!    M12   -    20    -    3.700
!    M13   -    22    -    4.050
!    M14   -    29    -    8.550
!    M15   -    31    -   10.763
!    M16   -    32    -   12.013
!
!--------------------------------------------------------------------------------------

module cx_read_ahi_mod
   
   use date_tools_mod 
   
   implicit none  
   private
   public :: get_ahi_data
   
   
   
   integer, parameter :: NUM_CHN = 16
   
   type, public :: ahi_config_type
      character ( len = 255 ) :: file_base
      character ( len = 400 ) :: data_path
      logical :: chan_on(16)
      character (len = 500) :: filename ( NUM_CHN)
      character ( len =20) :: varname (NUM_CHN)
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
      type ( ahi_chn_type ) :: chn (NUM_CHN)
      type ( geo_str ) :: geo
      type ( date_type ) :: time_start_obj
      type ( date_type ) :: time_end_obj
      contains
      procedure :: deallocate_all
      procedure :: allocate_geo
   end type ahi_data_out_type
   
          
contains
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine get_ahi_data ( config , out , only_nav )
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: out
      logical, optional, intent(in) :: only_nav
      
      call set_filenames ( config )
      call set_times ( config % file_base , out  )
      call read_navigation ( config , out )
      
      if ( .not. present ( only_nav )) then
         call read_ahi_level1b ( config , out )
      end if
      
   end subroutine get_ahi_data
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine set_times ( file_base , out )
      character ( len = * ) :: file_base 
      type ( ahi_data_out_type ) :: out
      
      integer :: year
      integer :: month
      integer :: day
      integer :: hour
      integer :: minute
      
      
      read(file_base(8:11), fmt="(I4)") year
      read(file_base(12:13), fmt="(I2)") month
      read(file_base(14:15), fmt="(I2)") day
      read(file_base(17:18), fmt="(I2)") hour
      read(file_base(19:20), fmt="(I2)") minute
      
      call out % time_start_obj % set_date ( &
            year , month, day , hour, minute &
            )
      
      out % time_end_obj = out % time_start_obj
      
      call out % time_end_obj % add_time ( minute = 8)
       
   end
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine set_filenames ( config)
      use string_functions, only: replace_text
      
      type ( ahi_config_type ) :: config
      
      integer :: i
      character(len=2) :: identifier
      character ( len=255) :: file_for_this_channel
      

      do i = 1 , 16
        
         write (identifier , fmt ='(i2.2)') i
         
         file_for_this_channel = replace_text ( config % file_base,'B01','B'//identifier)
         config % filename ( i ) = trim (config % data_path)//trim(file_for_this_channel) 
         config % varname ( i ) = '/RAD'
         
      end do

   end subroutine 

   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine read_navigation ( config, ahi )  
      use viewing_geometry_module
      use readh5dataset
      
      implicit none
      
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: ahi
      
      character (len=120) ::  attr_name
      integer :: status
      integer:: i_chn
      
      real (8) :: cfac
      real (8) :: coff
      real (8) :: lfac
      real (8) :: loff
      real (8) :: sub_lon 
      real (8) :: sub_lat  
      real (8) :: latx, lonx
      real (8) :: fargc
      integer  :: ii,jj
      integer  :: x_full_disk
      integer  :: y_full_disk
   
      real :: GEO_ALTITUDE = 35786.0 !km
      
      integer :: VALID_PROJECTION = 7
      
      integer :: day_of_year
      real :: hour_frac
      
       ! - navigation
      
      ! - should be removed later then we trust the values in the file! 
      cfac = 20466276.
      coff = 2750.
      lfac = 20466276.
      loff = 2750.
      sub_lon =  140.70
      sub_lat = 0.
      
      ! - read in the constants fom file

      attr_name = trim('Projection/CFAC')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION ) ) , trim ( attr_name ), CFAC )
      attr_name = trim('Projection/LFAC')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION) ) , trim ( attr_name ), LFAC )
      attr_name = trim('Projection/COFF')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION ) ) , trim ( attr_name ), COFF )
      attr_name = trim('Projection/LOFF')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION ) ) , trim ( attr_name ), LOFF )
      attr_name = trim('Projection/latitude_of_projection_origin')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION ) ) , trim ( attr_name ), sub_lat )
      attr_name = trim('Projection/longitude_of_projection_origin')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION) ) , trim ( attr_name ), sub_lon )

      call ahi % allocate_geo (config % h5_count(1),config % h5_count(2)) 
      
      call ahi % time_start_obj % get_date ( doy = day_of_year, hour_frac = hour_frac )

      do jj = 1 , config % h5_count(2)     
         do ii = 1 ,   config % h5_count(1)
            
            x_full_disk = ii + config % h5_offset(1)
            y_full_disk = jj + config % h5_offset(2)
            
            call fgf_to_earth ( 3, dble(x_full_disk), dble(y_full_disk) , cfac, coff, lfac, loff, sub_lon &
               , lonx , latx )
           
            ahi % geo % lat (ii,jj) = latx  
            ahi % geo % lon (ii,jj) = lonx  
            
            if ( lonx == -999. ) cycle
            
            call  possol ( day_of_year ,  hour_frac  , real(lonx) &
                  , real(latx),ahi % geo % solzen (ii,jj),ahi % geo % solaz (ii,jj) )
                        
            ahi % geo % satzen (ii,jj) = sensor_zenith ( GEO_ALTITUDE, real(sub_lon),real(sub_lat),real(lonx) , real(latx) )
             
            ahi % geo % sataz (ii,jj) = sensor_azimuth (  real(sub_lon),real(sub_lat),real(lonx) , real(latx) )
            
            ahi % geo % relaz(ii,jj) = relative_azimuth (ahi % geo % solaz (ii,jj) , ahi % geo % sataz (ii,jj))
           
            ahi % geo % glintzen(ii,jj) = glint_angle (ahi % geo % solzen (ii,jj) ,ahi % geo % satzen (ii,jj) &
                  ,  ahi % geo % relaz(ii,jj) )
           
            ahi % geo % scatangle(ii,jj) = scattering_angle ( ahi % geo % solzen (ii,jj)  &
                  , ahi % geo % satzen (ii,jj), ahi % geo % relaz(ii,jj))
  
         end do
      end do   
      
      where ( ahi % geo % lat == -999.0)
         ahi % geo % is_space = .false.
      end where

   end subroutine read_navigation
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine read_ahi_level1b ( config, ahi )
   
      use readh5dataset
      
      
      implicit none
      
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: ahi
     
      integer :: status
           
      integer(kind = 2), pointer :: i2d_buffer( : , : )
      integer:: i_chn
      character (len=120) :: attr_name
      real (8) :: scale_factor, add_offset
      integer ( kind = 2 ) :: fillvalue
      integer ( kind = 4 ) , allocatable :: buffer_fake_i4 (:,:)
      real ( 8 ) :: cprime 
      logical :: is_solar_channel = .false.
      real(8) :: fargc
      integer :: ii,jj
            
      ! - executable

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
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine allocate_geo ( this, nx , ny )
      class ( ahi_data_out_type ) :: this
      integer, intent(in) :: nx 
      integer, intent(in) :: ny
      
      call this % deallocate_all 
          
      allocate (  this % geo % lon        (nx , ny) )
      allocate (  this % geo % lat        (nx , ny) )
      allocate (  this % geo % solzen     (nx , ny) )
      allocate (  this % geo % solaz      (nx , ny) )
      allocate (  this % geo % satzen     (nx , ny) )
      allocate (  this % geo % sataz      (nx , ny) )
      allocate (  this % geo % relaz      (nx , ny) )
      allocate (  this % geo % glintzen   (nx , ny) )
      allocate (  this % geo % scatangle  (nx , ny) )
      allocate (  this % geo % is_space   (nx , ny) )
      
      this % geo % lon        =   -999.
      this % geo % lat        =   -999.
      this % geo % solzen     =   -999.
      this % geo % solaz      =   -999.
      this % geo % satzen     =   -999.
      this % geo % sataz      =   -999.
      this % geo % relaz      =   -999.
      this % geo % glintzen   =   -999.
      this % geo % scatangle  =   -999.
      this % geo % is_space = .true.
      
   
   end subroutine allocate_geo
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine deallocate_all (this )
      class ( ahi_data_out_type ) :: this
      integer :: i_chn
      if (allocated ( this % geo % lon)) deallocate ( this % geo % lon) 
      if (allocated ( this % geo % lat)) deallocate ( this % geo % lat) 
      
           
      if ( allocated  (  this % geo % solzen   ) ) deallocate (  this % geo % solzen   )
      if ( allocated  (  this % geo % solaz      ) ) deallocate (  this % geo % solaz   )
      if ( allocated  (  this % geo % satzen     ) ) deallocate (  this % geo % satzen  )
      if ( allocated  (  this % geo % sataz     ) ) deallocate (  this % geo % sataz  )
      if ( allocated  (  this % geo % relaz      ) ) deallocate (  this % geo % relaz   )
      if ( allocated  (  this % geo % glintzen    ) ) deallocate (  this % geo % glintzen   )
      if ( allocated  (  this % geo % scatangle   ) ) deallocate (  this % geo % scatangle   )
      if ( allocated  (  this % geo % is_space    ) ) deallocate (  this % geo % is_space  )
      
      do i_chn = 1 ,16
         if (allocated (  this  % chn(i_chn) % rad )) deallocate (this  % chn(i_chn) % rad)
         if (allocated (  this  % chn(i_chn) % bt ) ) deallocate (this  % chn(i_chn) % bt)
         if (allocated (  this  % chn(i_chn) % ref )) deallocate (this  % chn(i_chn) % ref)
      end do
   
   
   end subroutine deallocate_all

end module cx_read_ahi_mod
