! $Header$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: iff_module.f90 (src)
!       IFF_MODULE (program)
!
! PURPOSE: IFF read tool
!
! DESCRIPTION: This moduule reads IFF data
!
! AUTHORS:
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
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
!--------------------------------------------------------------------------------------
module IFF_MODULE

 use HDF
 use PLANCK
 use PIXEL_COMMON, only: &
     sensor &
     , image
 use CONSTANTS, only: &
     sym &
     , int1 &
     , int2 &
     , int4 &
     , int8 &
     , real4 & 
     , real8 & 
     , ipre
 use CALIBRATION_CONSTANTS, only: &
     Planck_Nu

 implicit none

 public :: GET_IFF_DATA
 public :: READ_IFF_DATE_TIME

   type , public :: iff_data_config
      integer , dimension( 2 ) :: year_int
      integer , dimension( 2 ) :: doy_int
      integer :: n_chan
      integer , dimension( 37 ) :: chan_list
      integer , dimension( 2 )  :: offset
      integer , dimension( 2 ) :: count
      logical , dimension( 37 ) :: chan_on
      logical :: iff_cloud_mask_on
      character ( len = 255 ) :: iff_file
      character ( len = 355 ) :: dir_1b
      character ( len = 355 ) :: Ancil_Data_Dir
   end type iff_data_config

   type :: geo_str
      real , dimension (:,:) , allocatable :: solzen
      real , dimension (:,:) , allocatable :: satzen
      real , dimension (:,:) , allocatable :: solaz
      real , dimension (:,:) , allocatable :: sataz
      real , dimension (:,:) , allocatable :: relaz
      real , dimension (:,:) , allocatable :: lat
      real , dimension (:,:) , allocatable :: lon
      real , dimension (:) , allocatable :: scan_time
      integer , dimension (:) , allocatable :: ascend
      integer(kind=int1) , dimension (:,:) , allocatable :: sounder_fov
      integer(kind=int2) , dimension (:,:) , allocatable :: sounder_x
      integer(kind=int2) , dimension (:,:) , allocatable :: sounder_y
   end type  geo_str

   type :: band_str
      logical :: is_read
      real,dimension (:,:) , allocatable :: ref
      real,dimension (:,:) , allocatable :: rad
      real,dimension (:,:) , allocatable :: bt
      contains
      procedure  ::  dealloc_band_str
   end type band_str

   type :: cloud_products_str
      integer , dimension(:,:) , allocatable :: cld_mask
      real , dimension(:,:) , allocatable :: sndr_cld_temp ! MJH
      real , dimension(:,:) , allocatable :: sndr_cld_pres
      real , dimension(:,:) , allocatable :: sndr_cld_height 
      integer , dimension(:,:) , allocatable :: sndr_mask 
      integer , dimension(:,:) , allocatable :: sndr_ele_idx
      integer , dimension(:,:) , allocatable :: sndr_line_idx 
   end type cloud_products_str

   type, public :: iff_data_out
      type (band_str), allocatable :: band(:)
      type (geo_str) :: geo
      type ( cloud_products_str) :: prd
      contains
      procedure  :: dealloc => dealloc_iff_data_out
   end type iff_data_out

   real(kind=real4), parameter, private:: missing_value = -999.0

 contains

!----------------------------------------------------------------------
!   Subroutine to get dimentions from IFF 1b file Latitude
!
subroutine GET_IFF_DIMS (File_Name, Nx, Ny)

   use HDF_READ_MODULE

   implicit none

   character(len=*), intent(in):: File_Name
   integer, intent(out):: Nx
   integer, intent(out):: Ny

   integer :: Status
   integer :: Id
   integer :: Rank
   integer, dimension (10) :: Sds_Dims
   character(len=64):: Sds_Name

   Status = 0
   Sds_Name = "Latitude"

   error_check: do while (Status == 0)

   Status = OPEN_FILE_HDF_READ ( TRIM(File_Name), Id ) + Status
   Status = HDF_SDS_DIMENIONS_READER ( Id, TRIM(Sds_Name), Rank, Sds_Dims ) + Status
   call CLOSE_FILE_HDF_READ( Id, TRIM(File_Name) )

   Nx = Sds_Dims(1)
   Ny = Sds_Dims(2)

   Status = 1
   enddo  error_check ! end of while loop

end subroutine GET_IFF_DIMS

!----------------------------------------------------------------------
!   Main subroutine to read IFF files
!
subroutine GET_IFF_DATA ( config , out )
      type ( iff_data_config ) , intent ( inout ) :: config
      type ( iff_data_out ) , intent ( out ) :: out
      integer :: error_out

      call CHECK_INPUT ( config )

      ! -- call
      call READ_IFF_LEVEL1B ( config, out , error_out )

end subroutine GET_IFF_DATA

!----------------------------------------------------------------------
!  check input - should be extended  ( testing if file exists etc...
!
   subroutine CHECK_INPUT ( config )                                                                                                       
       type ( iff_data_config ) , intent ( inout ) :: config

      if ( config % offset (1) <= 0 )  config % offset (1) = 1
      if ( config % offset (2) <= 0 )  config % offset (2)  = 1
      
      if ( config % count (1) <= 0 ) then
         select case ( trim(Sensor%Sensor_Name) )
         case ('VIIRS-IFF') 
            config % count (1) = 3200
         case ('AQUA-IFF') 
            config % count (1) = 1354
         case ( 'AVHRR-IFF') 
            config % count (1) = 409
         case default      
         end select   
      end if   
      if ( config % count (2) <= 0 )  config % count (2)  = Image%Number_Of_Lines_Per_Segment

   end subroutine CHECK_INPUT

!----------------------------------------------------------------------                                                                      
!   Subroutine to read IFF 1b files                                                                                                        
!
subroutine READ_IFF_LEVEL1B ( config, out, error_out )

      use FILE_TOOLS
      use HDF_READ_MODULE
      use NUMERICAL_ROUTINES, only: &
          COUNTSUBSTRING &
          , SPLIT_STRING &
          , REPLACE_CHAR_IN_STRG



      type ( iff_data_config ) , intent ( in ) :: config
      type ( iff_data_out ) , intent ( out ) :: out
      integer, intent(out) , optional :: error_out

      integer(kind=int4) :: iend
      integer(kind=int4) :: i_geo
      integer(kind=int4) :: i_band
      integer(kind=int4) :: iband_sds
      integer(kind=int4) :: ii_ref_rad
      integer(kind=int4) :: Status
      integer(kind=int4) :: Id
      integer(kind=int4) :: Rank
      integer(kind=int4) :: num_ref_ch
      integer(kind=int4) :: num_rad_ch
      integer(kind=int4) :: n_files
      integer(kind=int4) :: ii
      integer(kind=int4) :: num_char_band_names
      integer(kind=int4) :: nx_start , nx_end , ny_start , ny_end
      integer(kind=int4), dimension (10) :: Sds_Dims
      integer(kind=int4), dimension(:), allocatable  :: band_names_int_ref
      integer(kind=int4), dimension(:), allocatable  :: band_names_int_rad
      integer(kind=int4), dimension (2) :: start_2d, stride_2d, edge_2d
      integer(kind=int4), dimension (3) :: start_3d, stride_3d, edge_3d
      integer(kind=int1), dimension( : , : , : ) , allocatable :: i3d_buffer
      ! MJH Aug 2015 these are unnecessary now; use i2d_8_buffer instead
      integer(kind=int2), dimension( : , :) , allocatable :: i2d_buffer_int2
      integer(kind=int1), dimension( : , :) , allocatable :: i2d_buffer_int1
      integer(kind=int4) , dimension(2) :: dim_seg
      integer(kind=int1), dimension( : , : ) , allocatable :: i2d_8_buffer
      integer(kind=int2), dimension( : , : ) , allocatable :: i2d_16_buffer
      real(kind=int4), parameter :: fill_value = 65535.
      real(kind=int8), parameter :: sec_per_day = 86400.
      real(kind=int4), dimension(36) :: nu_list
      real(kind=int4), dimension ( : ) , allocatable:: time_msec_day
      real(kind=int8), dimension( : ) , allocatable :: r1d_buffer
      real(kind=int4), dimension( : , : ) , allocatable :: r2d_buffer
      real(kind=int4), dimension( : , : , : ) , allocatable :: r3d_buffer
      character (len=100), dimension ( 7 ) :: setname_geo_list = (/ character(len=40) :: &
                           'Latitude'            & ! 1
                         , 'Longitude'           & ! 2
                         , 'SensorAzimuth'       & ! 3
                         , 'SensorZenith'        & ! 4
                         , 'SolarAzimuth'        & ! 5
                         , 'SolarZenith'         & ! 6
                         , 'ScanStartTime'      /) ! 7
      character (len=100) :: sndr_cld_temp_strg = 'SounderCloudTemperature' !MJH
      character (len=100) :: sndr_cld_pres_strg = 'SounderCloudPressure' 
      character (len=100) :: sndr_cld_height_strg = 'SounderCloudHeights' 
      character (len=100) :: sndr_mask_strg = 'HirsMask' 
      character (len=100) :: sndr_ele_idx_strg = 'SounderX' 
      character (len=100) :: sndr_line_idx_strg = 'SounderY' 
      character (len = 150) :: setname_band
      character (len = 255) , pointer , dimension (:) :: file_arr_dummy
      character (len = 250) :: file_cld_mask
      character (len = 250) :: file_srch
      character (len = 200) :: band_names_char
      character (len = 100), dimension(:), allocatable :: band_names_char_arr
      character (len = 4) :: year_strg
      character (len = 3) :: doy_strg
 
      error_out = 0
      iend = 0
      Status = 0

      error_check: do while (Status == 0 .and. iend == 0)

      ! calculate start, end
      nx_start = config % offset(1) - 1
      ny_start = config % offset(2) - 1
      start_2d = [nx_start, ny_start]
      stride_2d = [ 1, 1 ]
      nx_end = nx_start + config % count(1)
      ny_end = ny_start + config % count(2)
      dim_seg = [ nx_end - nx_start, ny_end - ny_start ]
      edge_2d = dim_seg

      ! set nu numbes for radiance conversion
      nu_list = 0
      if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
         nu_list(20:25) = [Planck_Nu(20), Planck_Nu(20), Planck_Nu(31), Planck_Nu(23), Planck_Nu(24), Planck_Nu(25)]
         nu_list(27:36) = [Planck_Nu(27), Planck_Nu(28), Planck_Nu(32), Planck_Nu(30), Planck_Nu(31), Planck_Nu(32), &
                           Planck_Nu(33), Planck_Nu(34), Planck_Nu(35), Planck_Nu(36)]
      else
         nu_list(20:25) = [Planck_Nu(20), Planck_Nu(21), Planck_Nu(22), Planck_Nu(23), Planck_Nu(24), Planck_Nu(25)]
         nu_list(27:37) = [Planck_Nu(27), Planck_Nu(28), Planck_Nu(29), Planck_Nu(30), Planck_Nu(31), Planck_Nu(32), &
                           Planck_Nu(33), Planck_Nu(34), Planck_Nu(35), Planck_Nu(36), Planck_Nu(33)]
      endif

      ! open file
      Status = OPEN_FILE_HDF_READ ( TRIM( config % iff_file ), Id ) + Status
      Status = HDF_SDS_DIMENIONS_READER ( Id, TRIM(setname_geo_list(1)), Rank, Sds_Dims ) + Status

      ! --- read geo info
      do i_geo = 1 , 6
         Status = READ_HDF_SDS_FLOAT32_2D(Id,TRIM(setname_geo_list(i_geo)),start_2d, &
                     stride_2d,edge_2d,r2d_buffer) + Status

         select case (i_geo)
         case(1)
            if (.not. allocated ( out % geo % lat ) ) allocate ( out % geo % lat (dim_seg(1), dim_seg(2)) )
            out % geo % lat = r2d_buffer
         case(2)
            if (.not. allocated ( out % geo % lon ) ) allocate ( out % geo % lon (dim_seg(1), dim_seg(2)) )
            out % geo % lon = r2d_buffer
         case(3)
            if (.not. allocated ( out % geo % sataz ) ) allocate ( out % geo % sataz (dim_seg(1), dim_seg(2)) )
            out % geo % sataz = r2d_buffer
         case(4)
            if (.not. allocated ( out % geo % satzen ) ) allocate ( out % geo % satzen (dim_seg(1), dim_seg(2)) )
            out % geo % satzen = r2d_buffer
         case(5)
            if (.not. allocated ( out % geo % solaz ) ) allocate ( out % geo % solaz (dim_seg(1), dim_seg(2)) )
            out % geo % solaz = r2d_buffer
         case(6)
            if (.not. allocated ( out % geo % solzen ) ) allocate ( out % geo % solzen (dim_seg(1), dim_seg(2)) )
            out % geo % solzen = r2d_buffer
        end select

      end do ! read geo info

      ! MJH if AVHRR, read in aux fields like sounder collocation indices
      ! and Menzel HIRS cloud heights
      if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
          print *, 'Reading AVHRR-IFF aux fields'
          ! cld temp
          Status = READ_HDF_SDS_FLOAT32_2D(Id,TRIM(sndr_cld_temp_strg),start_2d, &
                      stride_2d,edge_2d,r2d_buffer) + Status
          if (.not. allocated ( out % prd % sndr_cld_temp ) ) allocate &
              ( out % prd % sndr_cld_temp (dim_seg(1), dim_seg(2)) )
          out % prd % sndr_cld_temp = r2d_buffer
          ! cld pressure
          Status = READ_HDF_SDS_FLOAT32_2D(Id,TRIM(sndr_cld_pres_strg),start_2d, &
                      stride_2d,edge_2d,r2d_buffer) + Status
          if (.not. allocated ( out % prd % sndr_cld_pres ) ) allocate &
              ( out % prd % sndr_cld_pres (dim_seg(1), dim_seg(2)) )
          out % prd % sndr_cld_pres = r2d_buffer
          ! cld height
          Status = READ_HDF_SDS_FLOAT32_2D(Id,TRIM(sndr_cld_height_strg),start_2d, &
                      stride_2d,edge_2d,r2d_buffer) + Status
          if (.not. allocated ( out % prd % sndr_cld_height ) ) allocate &
              ( out % prd % sndr_cld_height (dim_seg(1), dim_seg(2)) )
          out % prd % sndr_cld_height = r2d_buffer
          ! hirs mask
          Status = READ_HDF_SDS_INT8_2D(Id,TRIM(sndr_mask_strg),start_2d, &
                      stride_2d,edge_2d,i2d_buffer_int1) + Status
          if (.not. allocated ( out % prd % sndr_mask ) ) allocate &
              ( out % prd % sndr_mask (dim_seg(1), dim_seg(2)) )
          out % prd % sndr_mask = i2d_buffer_int1
          ! element index in original hirs granule
          Status = READ_HDF_SDS_INT16_2D(Id,TRIM(sndr_ele_idx_strg),start_2d, &
                      stride_2d,edge_2d,i2d_buffer_int2) + Status
          if (.not. allocated ( out % prd % sndr_ele_idx ) ) allocate &
              ( out % prd % sndr_ele_idx (dim_seg(1), dim_seg(2)) )
          out % prd % sndr_ele_idx = i2d_buffer_int2
          ! line index in original hirs granule
          Status = READ_HDF_SDS_INT16_2D(Id,TRIM(sndr_line_idx_strg),start_2d, &
                      stride_2d,edge_2d,i2d_buffer_int2) + Status
          if (.not. allocated ( out % prd % sndr_line_idx ) ) allocate &
              ( out % prd % sndr_line_idx (dim_seg(1), dim_seg(2)) )
          out % prd % sndr_line_idx = i2d_buffer_int2
      endif

      ! MJH Moved this out of geo loop 
      ! (it never gets allocated?? see https://software.intel.com/en-us/forums/topic/295774)
      deallocate ( r2d_buffer )
      deallocate ( i2d_buffer_int2 )
      deallocate ( i2d_buffer_int1 )

      ! read scan time and convert to milliseconds
      Status = READ_HDF_SDS_FLOAT64_1D(Id,TRIM(setname_geo_list(7)),start_2d(2), &
                     stride_2d(2),edge_2d(2),r1d_buffer) + Status
      allocate ( time_msec_day ( size ( r1d_buffer)))                                                                                               
      time_msec_day = ( dmod ( r1d_buffer , sec_per_day ) ) * 1000
      if (.not. allocated ( out % geo % scan_time ) ) allocate ( out % geo % scan_time (dim_seg(2)) )
      out % geo % scan_time = time_msec_day
      deallocate ( time_msec_day )
      deallocate ( r1d_buffer )

      ! --- Read reflectence band centers to find position in the array

      ! loop over 1 = reflectance and 2 = radiance
      do ii_ref_rad = 1, 2

         if ( ii_ref_rad == 1 ) then
            setname_band = 'ReflectiveSolarBands'
         elseif ( ii_ref_rad == 2 ) then
            setname_band = 'EmissiveBands'
         endif

         ! get band names from the attribute
         Status = READ_HDF_ATTRIBUTE_CHAR8_SCALAR(Id, TRIM(setname_band), 'band_names', band_names_char) + Status

         ! find out how many attributes
         num_char_band_names = COUNTSUBSTRING ( band_names_char, ',' ) + 1
         ! split one string to channel names array
         Status = SPLIT_STRING (band_names_char, ',', num_char_band_names, band_names_char_arr) + Status

         ! delete all unuseful characters from the string
         ! for CrIS, HIRS and AIRS channels add 1 (131 - 136)
         ! for MODIS set low 13 & 14 channels to 913 and 914 to ignore it
         do ii = 1, num_char_band_names
            Status = REPLACE_CHAR_IN_STRG (band_names_char_arr(ii),'-','1','before') + Status ! Sounder
            Status = REPLACE_CHAR_IN_STRG (band_names_char_arr(ii),'o','9','before') + Status ! Pseudo
            select case (trim(Sensor%Sensor_Name))
            case('VIIRS-IFF')
               Status = REPLACE_CHAR_IN_STRG (band_names_char_arr(ii),'M',' ','before') + Status ! M-Bands
            case('AQUA-IFF')
               Status = REPLACE_CHAR_IN_STRG (band_names_char_arr(ii),'l','9','after') + Status ! lo/hi
               Status = REPLACE_CHAR_IN_STRG (band_names_char_arr(ii),'h',' ','after') + Status  
            case('AVHRR-IFF')
               Status = REPLACE_CHAR_IN_STRG (band_names_char_arr(ii),'a',' ','after') + Status ! 3a/3b
               Status = REPLACE_CHAR_IN_STRG (band_names_char_arr(ii),'b',' ','after') + Status               
            case default
            
            end select  
         
         enddo ! loop over band names

         ! convert string array to integer
         if ( ii_ref_rad == 1 ) then
            allocate ( band_names_int_ref (num_char_band_names) )
            read (band_names_char_arr,fmt=*) band_names_int_ref
            num_ref_ch = num_char_band_names
            if (trim(Sensor%Sensor_Name) == 'VIIRS-IFF') then
               do ii = 1, num_char_band_names
                 select case ( band_names_int_ref(ii) )
                    case (1)
                       band_names_int_ref(ii) = 8
                    case (2)
                       band_names_int_ref(ii) = 9
                    case (3)
                       band_names_int_ref(ii) = 3
                    case (4)
                       band_names_int_ref(ii) = 4
                    case (5)
                       band_names_int_ref(ii) = 1
                    case (6)
                       band_names_int_ref(ii) = 15
                    case (7)
                       band_names_int_ref(ii) = 2
                    case (8)
                       band_names_int_ref(ii) = 5
                    case (9)
                       band_names_int_ref(ii) = 26
                    case (10)
                       band_names_int_ref(ii) = 6
                    case (11)
                       band_names_int_ref(ii) = 7
                 end select
               enddo
            elseif (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
               do ii = 1, num_char_band_names
                 select case ( band_names_int_ref(ii) )
                    case (1)
                       band_names_int_ref(ii) = 1
                    case (2)
                       band_names_int_ref(ii) = 2
                    case (3)
                       band_names_int_ref(ii) = 6
                 end select
               enddo
            endif
         elseif ( ii_ref_rad == 2 ) then
            allocate ( band_names_int_rad (num_char_band_names) )
            read (band_names_char_arr,fmt=*) band_names_int_rad
            num_rad_ch = num_char_band_names
            if (trim(Sensor%Sensor_Name) == 'VIIRS-IFF') then
               do ii = 1, num_char_band_names
                  select case ( band_names_int_rad(ii) )
                     case (12) ! M-Bands
                        band_names_int_rad(ii) = 20
                     case (13) ! M-Bands
                        band_names_int_rad(ii) = 22
                     case (14) ! M-Bands
                        band_names_int_rad(ii) = 29
                     case (15) ! M-Bands
                        band_names_int_rad(ii) = 31
                     case (16) ! M-Bands
                        band_names_int_rad(ii) = 32
                     case (133) ! Sounder
                        band_names_int_rad(ii) = 33
                     case (134) ! Sounder
                        band_names_int_rad(ii) = 34
                     case (135) ! Sounder
                        band_names_int_rad(ii) = 35
                     case (136) ! Sounder
                        band_names_int_rad(ii) = 36
                     case (933) ! Pseudo
                        band_names_int_rad(ii) = 45
                  end select
               enddo
            elseif (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
               do ii = 1, num_char_band_names
                  select case ( band_names_int_rad(ii) )
                     case (3)
                        band_names_int_rad(ii) = 20
                     case (4)
                        band_names_int_rad(ii) = 31
                     case (5)
                        band_names_int_rad(ii) = 32
                     case (14) ! HIRS-4
                        band_names_int_rad(ii) = 36
                     case (15) ! HIRS-5
                        band_names_int_rad(ii) = 35
                     case (16) ! HIRS-6
                        band_names_int_rad(ii) = 34
                     case (17) ! HIRS-7
                        band_names_int_rad(ii) = 33
                     ! Attn: Use unused ch. to save HIRS 11um
                     case (18) ! HIRS-8 ! same as AVHRR ch4
                        band_names_int_rad(ii) = 22
                     case (19) ! HIRS-9
                        band_names_int_rad(ii) = 30
                     ! Attn: Use ch. 29 to save either HIRS 12um or 8.55um
                     case (110) ! HIRS-10  ! same as AVHRR ch5
                        band_names_int_rad(ii) = 29 ! early 29, latter 32
                     case (111) ! HIRS-11
                        band_names_int_rad(ii) = 28
                     case (112) ! HIRS-12
                        band_names_int_rad(ii) = 27
                     case (114) ! HIRS-14
                        band_names_int_rad(ii) = 25
                     case (116) ! HIRS-16
                        band_names_int_rad(ii) = 24
                     case (118) ! HIRS-18
                        band_names_int_rad(ii) = 23
                     ! Attn: Use unused ch. to save HIRS 3.75um
                     case (119) ! HIRS-19 ! same as AVHRR ch3
                        band_names_int_rad(ii) = 21
                     case (933) ! Pseudo
                        band_names_int_rad(ii) = 45
                  end select
               enddo
            elseif (trim(Sensor%Sensor_Name) == 'AQUA-IFF') then
               do ii = 1, num_char_band_names                                                                                                                                      
                  if ( band_names_int_rad(ii) == 933 ) band_names_int_rad(ii) = 45
               enddo
            endif
         endif

         ! --- dealocate variables
         deallocate ( band_names_char_arr )

      ! --- end loop over ref & rad
      enddo
      
      allocate ( out % band (config % n_chan))
      ! --- read bands
      do i_band = 1, config % n_chan

         ! check if channel is on
         if ( .not. config % chan_on ( i_band ) ) cycle

         iband_sds = - 1
         ! based on channel, set appropriate names and sds band index
         if ((config % chan_list(i_band) >= 1 .and. config % chan_list(i_band) <= 19) &
            .or. config % chan_list(i_band) == 26) then
            setname_band = 'ReflectiveSolarBands'

            ! find what current channel has array number
            do ii = 1, num_ref_ch
               if ( config % chan_list(i_band) == band_names_int_ref(ii) ) then
                  iband_sds = ii
               endif
            enddo

         elseif ((config % chan_list(i_band) >= 20 .and. config % chan_list(i_band) <= 36 &
            .and. config % chan_list(i_band) /= 26) .or. config % chan_list(i_band) == 45) then
            setname_band = 'EmissiveBands'

            ! find what current channel has array number
            do ii = 1, num_rad_ch
               if ( config % chan_list(i_band) == band_names_int_rad(ii) ) then
                  iband_sds = ii
               endif
            enddo

         endif

         ! read bands
         stride_3d = (/ 1, 1, 1 /)
         start_3d = (/ 0, ny_start, iband_sds - 1 /)
         edge_3d = (/ dim_seg(1), dim_seg(2), 1 /)
         allocate ( r3d_buffer(dim_seg(1), dim_seg(2), 1) )

         Status = READ_HDF_SDS_FLOAT32_3D(Id,TRIM(setname_band),start_3d, &
                     stride_3d,edge_3d,r3d_buffer) + Status

         ! change fill and unrealistic values to missing
         where(r3d_buffer == fill_value .or. r3d_buffer .gt. 400.)
            r3d_buffer = missing_value
         endwhere

         ! --- calibrate reflective channel data from 0 to 100%
         ! --- convert radiance and save all data
         if (trim(setname_band) == 'ReflectiveSolarBands') then
            r3d_buffer =  100.0 * r3d_buffer
            where(r3d_buffer .LT. 0.)
               r3d_buffer = missing_value
            endwhere
            if (.not. allocated ( out % band (i_band) % ref ) ) &
                      allocate (out % band (i_band) % ref (dim_seg(1), dim_seg(2)) )
            out % band (i_band) % ref =  r3d_buffer(:,:,1)
         elseif (trim(setname_band) == 'EmissiveBands') then
            call CONVERT_RADIANCE(r3d_buffer(:,:,1),nu_list(i_band),missing_value)
            if (.not. allocated ( out % band (i_band) % rad ) ) &
                     allocate (out % band (i_band) % rad (dim_seg(1), dim_seg(2)) )
            out % band (i_band) % rad =  r3d_buffer(:,:,1)

         endif
            
         out % band (i_band) % is_read = .true.

         ! --- dealocate variables
         deallocate ( r3d_buffer )

      enddo ! loop over channals

      deallocate ( band_names_int_ref )                                                                                        
      deallocate ( band_names_int_rad )


      ! --- Read indices from sounder
      setname_band = 'SounderFOV'
      stride_2d = (/ 1, 1 /)
      start_2d = (/ 0, ny_start /)
      edge_2d = (/ dim_seg(1), dim_seg(2) /)
      allocate ( i2d_8_buffer(dim_seg(1), dim_seg(2)) )
      allocate ( i2d_16_buffer(dim_seg(1), dim_seg(2)) )
      allocate ( out % geo % sounder_fov(dim_seg(1), dim_seg(2)) )
      allocate ( out % geo % sounder_x(dim_seg(1), dim_seg(2)) )
      allocate ( out % geo % sounder_y(dim_seg(1), dim_seg(2)) )

      Status = READ_HDF_SDS_INT8_2D(Id,TRIM(setname_band),start_2d, &
                  stride_2d,edge_2d,i2d_8_buffer) + Status
      out % geo % sounder_fov = i2d_8_buffer

      ! change fill values to missing
      where(i2d_8_buffer .lt. 0)
         out % geo % sounder_fov = -128
      endwhere

      setname_band = 'SounderX'
      Status = READ_HDF_SDS_INT16_2D(Id,TRIM(setname_band),start_2d, &
                  stride_2d,edge_2d,i2d_16_buffer) + Status
      out % geo % sounder_x = i2d_16_buffer

      ! change fill values to missing
      where(i2d_16_buffer .lt. 0)
         out % geo % sounder_x = -32768
      endwhere

      i2d_16_buffer = -1
      setname_band = 'SounderY'
      Status = READ_HDF_SDS_INT16_2D(Id,TRIM(setname_band),start_2d, &
                  stride_2d,edge_2d,i2d_16_buffer) + Status
      out % geo % sounder_y = i2d_16_buffer

      ! change fill values to missing
      where(i2d_16_buffer .lt. 0)
         out % geo % sounder_y = -32768
      endwhere

      ! --- dealocate variables
      deallocate ( i2d_8_buffer )
      deallocate ( i2d_16_buffer )

      !Close main file
      call CLOSE_FILE_HDF_READ(Id,TRIM(config%iff_file))

      !0        1         2         3         4         5         6
      !123456789012345678901234567890123456789012345678901234567890
      !IFFSDR_aqua_d20130113_t235500_c20150202152558_ssec_dev.hdf
      !MYD35_L2.A2013013.2355_aqua_IFF_v2_2_6.hdf
      !IFFSDR_npp_d20130113_t235500_c20150201214538_ssec_dev.hdf
      !MYD35_L2.A2013013.2355_npp_IFF_v2_2_6.hdf

      ! --- Read cld_mask_aux if it set in default file
      if ( config %  iff_cloud_mask_on ) then
         if (ny_end == dim_seg(2)) &

          print *, EXE_PROMPT, "Reading in AUX cloud mask"

          ! ---
          write (year_strg, "(I4)") config % year_int
          write (doy_strg, "(I0.3)") config % doy_int

          if (trim(Sensor%Sensor_Name)== 'VIIRS-IFF') then
             file_srch = 'MYD35_L2.A'//trim(year_strg)//trim(doy_strg)//'.'// &
                         trim(config % iff_file &
                         (len(trim(config % dir_1b))+23:len(trim(config % dir_1b))+26)) &
                         //'_npp_IFF*.hdf'
          endif
          if (trim(Sensor%Sensor_Name)== 'AQUA-IFF' ) then
             file_srch = 'MYD35_L2.A'//trim(year_strg)//trim(doy_strg)//'.'// &
                         trim(config % iff_file &
                         (len(trim(config % dir_1b))+24:len(trim(config % dir_1b))+27)) &
                         //'_aqua_IFF*.hdf'
          endif

         file_arr_dummy => FILE_SEARCH ( trim(config %dir_1b), file_srch, n_files )


         ! if file found read cloud mask
         if ( n_files /= 0 ) then
            file_cld_mask = file_arr_dummy(1)
            setname_band = 'Cloud_Mask'

            allocate ( i3d_buffer ( dim_seg(1), dim_seg(2), 1 ) )
            stride_3d = (/ 1, 1, 1 /)
            start_3d = (/ 0, ny_start, 0 /)
            edge_3d = (/ dim_seg(1), dim_seg(2), 1 /)

            ! open cloud mask file
            Status = OPEN_FILE_HDF_READ ( TRIM(config % dir_1b) &
                            // TRIM( file_arr_dummy(1) ), Id ) + Status

            ! read data
            Status = READ_HDF_SDS_INT8_3d(Id,TRIM(setname_band),start_3d, &
                     stride_3d,edge_3d,i3d_buffer) + Status

            ! extract cloud mask bits (# 1 & 2 of the first byte)
            if (.not. allocated ( out % prd % cld_mask ) ) &
                      allocate ( out % prd % cld_mask ( dim_seg(1), dim_seg(2) ) )
            out % prd % cld_mask = 0

            where ( btest ( i3d_buffer ( : , : , 1 ), 1 ) )
               out % prd % cld_mask  =   out % prd % cld_mask  + 1
            end where
            where ( btest ( i3d_buffer ( : , : , 1 ), 2 ) )
               out % prd % cld_mask  =   out % prd % cld_mask  + 2
            end where

            ! convert to clavr-x system (0 - confidently clear,
            ! 1 - prob. clear, 2 - prob. cloudy, 3 - conf. cloudy)
            out % prd % cld_mask = 3 - out % prd % cld_mask

            ! set pixels where cloud mask wasn't determind to missing
            where ( btest ( i3d_buffer ( : , : , 1 ), 0 ) .eqv. .false. )
               out % prd % cld_mask  = missing_value 
            end where

            deallocate ( i3d_buffer)

            !Close cloud mask file
            call CLOSE_FILE_HDF_READ( Id, TRIM(config %dir_1b) &
                                    // TRIM(file_arr_dummy(1)) )

         else
            print*,'AUX Mask file not found , cld_mask_aux is set to missing'

         end if 

       endif

      iend = 1

      enddo  error_check ! end of while loop

end subroutine READ_IFF_LEVEL1B

!---------------------------------------------------------------------- 
!   Get the date & time from the file's name
!

subroutine READ_IFF_DATE_TIME(Path,Infile,year,doy,start_time, &
                  end_year,end_doy,end_time)
   
 use NUMERICAL_ROUTINES, only: &
     LEAP_YEAR_FCT &
     , JULIAN
   implicit none
   CHARACTER(Len=*), INTENT(IN) :: Path
   CHARACTER(Len=*), INTENT(IN) :: Infile
   INTEGER, INTENT(OUT) :: year
   INTEGER, INTENT(OUT) :: doy    !day of year
   INTEGER, INTENT(OUT) :: start_time  !millisec
   INTEGER, INTENT(OUT) :: end_time    !millisec
   INTEGER, INTENT(OUT) :: end_year
   INTEGER, INTENT(OUT) :: end_doy    !day of year

   integer(kind=int4):: Status_Flag
   integer(kind=int4):: Sd_Id
   integer(kind=int4):: Sds_Id
   integer(kind=int4):: DFACC_read
   integer(kind=int4):: sfend
   integer(kind=int4):: sfstart
   integer(kind=int4):: sfrcatt
   integer(kind=int4):: sffattr
   character(len=20):: Metadata_Temp

   INTEGER :: ileap
   INTEGER :: month
   INTEGER :: day
   INTEGER :: start_hour
   INTEGER :: start_minute
   INTEGER :: start_sec
   INTEGER :: end_hour
   INTEGER :: end_minute
   INTEGER :: end_sec
   INTEGER :: days_in_year

!print*,'start'
Status_Flag = 0

!---- open file
Sd_Id = sfstart(trim(Path)//trim(Infile), DFACC_read)
!print*,sd_id,trim(Path)//trim(Infile), DFACC_read
!--- determine attribute index                                                                                                 
Sds_Id = sffattr(Sd_Id, "granule_start_time")

!--- read start attribute
Status_Flag = sfrcatt(Sd_Id, Sds_Id, Metadata_Temp) + Status_Flag
!print*,sd_id,status_flag
!0        1         2
!12345678901234567890
!2013-01-13 00:00:00
!2013-01-13 00:04:59

if ( status_flag /= 0) then
   print*,sd_id,status_flag
   print*,'IFF MODULE: error reading granule'
   stop
end if
! --- Read data 
read(Metadata_Temp(1:4), fmt="(I4)") year
read(Metadata_Temp(6:7), fmt="(I2)") month
read(Metadata_Temp(9:10), fmt="(I2)") day
read(Metadata_Temp(12:13), fmt="(I2)") start_hour
read(Metadata_Temp(15:16), fmt="(I2)") start_minute
read(Metadata_Temp(18:19), fmt="(I2)") start_sec

! --- Calculate the date of year
ileap = 0
ileap = leap_year_fct(year)
days_in_year = 365 + ileap

!--- compute day of the year based on month
call JULIAN(day,month,year,doy)

! --- read end time
!--- determine attribute index
Sds_Id = sffattr(Sd_Id, "granule_end_time")

!--- read end attribute
Status_Flag = sfrcatt(Sd_Id, Sds_Id, Metadata_Temp) + Status_Flag

read(Metadata_Temp(12:13), fmt="(I2)") end_hour
read(Metadata_Temp(15:16), fmt="(I2)") end_minute
read(Metadata_Temp(18:19), fmt="(I2)") end_sec

!--- close file
Status_Flag = sfend(Sd_Id) + Status_Flag

! --- Calculate start END time
start_time = ((start_hour * 60 + start_minute) * 60 + start_sec) * 1000
end_time = ((end_hour * 60 + end_minute) * 60 + end_sec) * 1000

! --- Check if end time is in the next day or year
end_year = year
end_doy = doy
if ( end_time <= start_time) then
   end_doy = end_doy + 1
   if (end_doy > days_in_year) then
      end_year = end_year + 1
      end_doy = 1
   endif
endif

end subroutine READ_IFF_DATE_TIME

!----------------------------------------------------------------------
!  public routine to deallocate output structure
!  history: 11/14/2013 Denis B.
! 
!----------------------------------------------------------------------
   subroutine dealloc_iff_data_out ( this )
      class ( iff_data_out) :: this
      integer :: m

      if (allocated ( this%prd%cld_mask )) deallocate ( this%prd%cld_mask )
      if (allocated ( this%prd%sndr_cld_temp )) deallocate ( this%prd%sndr_cld_temp ) ! MJH
      if (allocated ( this%prd%sndr_cld_pres )) deallocate ( this%prd%sndr_cld_pres )
      if (allocated ( this%prd%sndr_cld_height )) deallocate ( this%prd%sndr_cld_height )
      if (allocated ( this%prd%sndr_mask )) deallocate ( this%prd%sndr_mask )
      if (allocated ( this%prd%sndr_ele_idx )) deallocate ( this%prd%sndr_ele_idx )
      if (allocated ( this%prd%sndr_line_idx )) deallocate ( this%prd%sndr_line_idx )

      if (allocated ( this%geo%solzen )) deallocate ( this%geo%solzen )
      if (allocated ( this%geo%satzen )) deallocate ( this%geo%satzen )
      if (allocated ( this%geo%solaz )) deallocate ( this%geo%solaz )
      if (allocated ( this%geo%sataz )) deallocate ( this%geo%sataz )
      if (allocated ( this%geo%relaz )) deallocate ( this%geo%relaz )
      if (allocated ( this%geo%lat )) deallocate ( this%geo%lat )
      if (allocated ( this%geo%lon )) deallocate ( this%geo%lon )
      if (allocated ( this%geo%scan_time )) deallocate ( this%geo%scan_time )
      if (allocated ( this%geo%sounder_fov )) deallocate ( this%geo%sounder_fov )
      if (allocated ( this%geo%sounder_x )) deallocate ( this%geo%sounder_x )
      if (allocated ( this%geo%sounder_y )) deallocate ( this%geo%sounder_y )

      do m = 1 , 37
        if (allocated (this%band (m) %ref )) deallocate ( this%band (m) %ref )
        if (allocated (this%band (m) %rad )) deallocate ( this%band (m) %rad )
        if (allocated (this%band (m) %bt )) deallocate ( this%band (m) %bt )
      end do
      
      if ( allocated ( this % band )) deallocate ( this % band)

   end subroutine dealloc_iff_data_out

!----------------------------------------------------------------------

   subroutine dealloc_band_str ( this )
      class ( band_str ) :: this
      print*,'band deallocation'
      if ( allocated (this%ref) ) deallocate (this%ref )
      if ( allocated (this%rad) ) deallocate (this%rad )
      if ( allocated (this%bt) ) deallocate (this%bt )
      print*,'done'
   end subroutine dealloc_band_str

!----------------------------------------------------------------------

end module IFF_MODULE

