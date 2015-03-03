!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: viirs_read_mod.f90 (src)
!       viirs_read_mod (program)
!
! PURPOSE: VIIRS read tool
!
! DESCRIPTION:  This module deals with reading VIIRS data
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
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
! HISTORY:   created      March 2013 (AW)
!                         avoids reading in full viirs data set, which caused seg 
!                         faults on some machines with low memory space ( 17 Oct 2013 AW ) 
!--------------------------------------------------------------------------------------
module viirs_read_mod

   use hdf5
   use readh5dataset
   implicit none
   
   private
   public :: GET_VIIRS_DATA
   public :: READ_NUMBER_OF_SCANS_FROM_VIIRS
   public :: READ_VIIRS_DATE_TIME_ATT
  
   integer, parameter, public:: int1 = selected_int_kind(1)
   integer, parameter, public:: int2 = selected_int_kind(3)
   integer, parameter, public:: int4 = selected_int_kind(8)
   integer, parameter, public:: int8 = selected_int_kind(10)
   integer, parameter, public:: real4 = selected_real_kind(6,37)
   integer, parameter, public:: real8 = selected_real_kind(15,307)
   integer, parameter, public:: ipre = real4
   
   ! - bowtie gaps values
   integer, parameter :: Ny_Pattern = 48
   integer, parameter :: Nx_Pattern = 3200
   logical, dimension(:,:), allocatable, public:: Gap_Pixel_Mask_pattern
   integer(kind=int4), dimension(:,:), allocatable, public:: Gap_Line_Idx_pattern
    
   type :: syms
      integer::  NO  =  0
      integer :: YES = 1
   end type syms
   type(syms),save ::sym
  
   logical :: gap_pattern_computed =.false.
   
   type , public :: viirs_data_config
      logical :: chan_on_rfl_mband ( 16 )
      logical :: chan_on_rad_mband ( 16 )
      logical :: chan_on_iband ( 5 )
      logical :: chan_on_dnb 
      logical :: viirs_cloud_mask_on
      logical :: viirs_cloud_type_on
      character ( len = 255 ) :: gmtco_file
      integer , dimension( 2 )  :: offset
      integer , dimension( 2 ) :: count 
      character ( len =355 ) :: dir_1b
      character ( len =255 ) :: file_gmtco_base
      character (len = 355 ) :: Ancil_Data_Dir
      real, dimension(16) :: Nu_List
   end type viirs_data_config
   
   type :: geo_str
      real , dimension (:,:) , allocatable :: solzen
      real , dimension (:,:) , allocatable :: satzen
      real , dimension (:,:) , allocatable :: solaz
      real , dimension (:,:) , allocatable :: sataz  
      real , dimension (:,:) , allocatable :: relaz 
      real , dimension (:,:) , allocatable :: lunaz 
      real , dimension (:,:) , allocatable :: lunzen 
      real , dimension (:,:) , allocatable :: lunrelaz
      real , dimension (:,:) , allocatable :: lat
      real , dimension (:,:) , allocatable :: lon  
      real , dimension (:) , allocatable :: scan_time
      integer , dimension (:) , allocatable :: ascend
      real :: Moon_Illum_Frac
      real :: Moon_Phase_Angle
   end type  geo_str
   
   type :: mband_str
      logical :: is_read
      real,dimension (:,:) , allocatable :: ref
      real,dimension (:,:) , allocatable :: rad
      real,dimension (:,:) , allocatable :: bt
      contains
      procedure  ::  dealloc_mband_str
   end type mband_str
   
   type :: dnb_on_mband_str
      real,dimension (:,:) , allocatable :: ref
      real,dimension (:,:) , allocatable :: rad
   end type dnb_on_mband_str

   type :: iband_str
      logical :: is_read
      real,dimension (:,:) , allocatable :: ref
      real,dimension (:,:) , allocatable :: bt
   end type iband_str
   
   type :: dnb_str
      real,dimension (:,:) , allocatable :: ref_lun
      real,dimension (:,:) , allocatable :: ref_sol
      real,dimension (:,:) , allocatable :: rad
   end type dnb_str
   
   type :: cloud_products_str
      integer , dimension(:,:) , allocatable :: cld_mask
      integer , dimension(:,:) , allocatable :: cld_phase
      integer , dimension(:,:) , allocatable :: cld_type
   end type cloud_products_str
   
   type :: gap_str
      logical , dimension ( : , : ) , allocatable :: mask
       
   end type gap_str

   type :: viirs_file_exists
      logical :: svm_file_exists ( 16 )
      logical :: svi_file_exists ( 5 )
      logical :: iicmo_file_exists
      logical :: svdnb_file_exists
      logical :: gdnbo_file_exists
   end type
  
   type, public :: viirs_data_out
      type (mband_str), dimension(16) :: mband
      type  (iband_str), dimension(5) :: iband    
      type ( geo_str) :: geo
      type (dnb_str) :: dnb
      type ( dnb_on_mband_str ) :: dnb_mgrid
      type ( cloud_products_str) :: prd
      type ( gap_str ) :: gap
      type ( viirs_file_exists ) :: file_exists
     
      contains
      procedure  :: dealloc =>dealloc_viirs_data_out
   end type viirs_data_out
   
contains
   !
   ! entree point routine
   !  
   subroutine get_viirs_data ( config , out )
      type ( viirs_data_config ) , intent ( inout ) :: config
      type ( viirs_data_out ) , intent ( out ) :: out
      integer :: error_out
      
      call check_input (config )
     
      ! -- call
      call read_viirs_level1b ( config, out , error_out)
     
   
   end subroutine get_viirs_data
   
   
   !
   !  check input - should be extended  ( testing if file exists etc...
   !
   subroutine check_input ( config )
       type ( viirs_data_config ) , intent ( inout ) :: config
   
      if ( config % offset (1) <= 0 )  config % offset (1) = 1
      if ( config % offset (2) <= 0 )  config % offset (2)  = 1
      if ( config % count (1) <= 0 )  config % count (1) = 3200
      if ( config % count (2) <= 0 )  config % count (2)  = 200
   
   end subroutine check_input
   
   

   !
   !
   !
   subroutine read_viirs_level1b (config, out,  error_out )
      use hdf5
      use readh5dataset
      use file_tools
   
      type ( viirs_data_config ) , intent ( in ) :: config
      type ( viirs_data_out ) , intent ( out ) :: out
      integer, intent(out) , optional :: error_out
           
      character(len=150) :: orbit_identifier
      integer(kind=int4) :: nx_start , nx_end , ny_start , ny_end
      integer(kind=int4) :: nx_start_iband  , ny_start_iband 
   
      real, dimension( : , : ) , pointer :: r2d_buffer
      real, dimension( : , : ) , allocatable :: r2d_data
      integer ( kind = int8) , dimension(:) , pointer :: i1d_buffer
      integer, dimension ( : , : ) , pointer :: i2d_buffer             
      integer, dimension ( :, :), allocatable :: cld_type_idps
      character ( len = 255 ) , pointer , dimension ( :) :: file_arr_dummy 
      
      character (len=100), dimension ( 7 ) :: setname_gmtco_list = (/ character (len =300) :: &
                          'All_Data/VIIRS-MOD-GEO-TC_All/Latitude              ' & ! 1
                         , 'All_Data/VIIRS-MOD-GEO-TC_All/Longitude            ' & ! 2
                         , 'All_Data/VIIRS-MOD-GEO-TC_All/StartTime            ' & ! 3
                         , 'All_Data/VIIRS-MOD-GEO-TC_All/SatelliteAzimuthAngle' & ! 4
                         , 'All_Data/VIIRS-MOD-GEO-TC_All/SatelliteZenithAngle ' & ! 5
                         , 'All_Data/VIIRS-MOD-GEO-TC_All/SolarAzimuthAngle    ' & ! 6
                         , 'All_Data/VIIRS-MOD-GEO-TC_All/SolarZenithAngle     ' & ! 7
                                                    /)                            
      integer :: i_gmtco
      character ( len =650 ) :: file_gmtco
      
      integer :: i_mband
      character ( len =100) :: setname_mband , setname_mband_fac
      character ( len =150 ) :: file_mband
      logical, dimension ( 16) :: is_mband_on
   
      integer :: i_iband
      character ( len =100) :: setname_iband, setname_iband_fac
      character ( len =150 ) :: file_iband
      logical, dimension ( 5) :: is_iband_on
   
      character ( len =150 ) :: file_dnb_idx
      character ( len =150 ) :: file_gdnbo, file_svdnb
      integer :: io_err_stat
      integer :: i_map
      integer , dimension (3200) :: d2m_indx
      
      character ( len =100) :: setname_iicmo
      character ( len =150 ) :: file_iicmo
      integer :: lun
      integer :: n_files
      integer :: k
      
      character ( len = 2 ) :: band_nr_file
      character ( len = 2 ) :: band_nr_var
      real , dimension(:) , pointer :: factors
      logical , dimension(16) :: data_scaled_mband
     
      integer, dimension(2) :: dim_buf
      logical, dimension(:,:), allocatable :: invalid_pixel
      integer(kind=int4) , dimension(2) :: dim_seg 
      integer(kind=int4) , dimension(2) :: dim_seg_iband
      integer(kind=int4) , dimension(2) :: dim_seg_dnb
      ! - this is to delete
      integer :: i
      integer , dimension (: ) , allocatable:: time_msec_day
      integer (kind = int8) , parameter :: microsec_per_day =  86400000000
      integer, dimension(2)::shape_buffer
      
      real :: missing_value_real4 = -999.
      integer :: missing_value_int = -999
     
      integer:: scaled_Missing 
      integer :: unscaled_Missing 
      integer , dimension(2) :: offset_mband , offset_iband
      integer :: day_of_year
      
      
      ! -- executable
    
      Scaled_Missing = 65528
      error_out = 0
      Unscaled_Missing = -999.3
 
      nx_start = config % offset(1)
      ny_start = config % offset(2)
      nx_end = nx_start + config % count(1) - 1
      ny_end = ny_start + config % count(2) - 1
      dim_seg = [ nx_end - nx_start +1, ny_end - ny_start + 1 ]
      
            
      file_gmtco = trim(config % dir_1b) // trim(config % file_gmtco_base)
      offset_mband = [ nx_start -1  , ny_start - 1] 
      do i_gmtco = 1 , 7 
        
         if ( i_gmtco == 3 ) then
             call h5readdataset ( file_gmtco, trim(setname_gmtco_list(i_gmtco)), i1d_buffer )
         else 
            call h5readdataset ( file_gmtco, trim(setname_gmtco_list(i_gmtco)), offset_mband, dim_seg, r2d_buffer )
            where ( r2d_buffer < unscaled_missing ) r2d_buffer = missing_value_real4
         end if
            
         select case (i_gmtco)
         case(1) 
            if (.not. allocated ( out % geo % lat) ) allocate ( out % geo % lat (dim_seg(1), dim_seg(2)) )
            out % geo % lat = r2d_buffer 
         case(2) 
            if (.not. allocated ( out % geo % lon) ) allocate ( out % geo % lon (dim_seg(1), dim_seg(2)) )
            out % geo %lon = r2d_buffer 
         case(3) 
            allocate ( time_msec_day ( size ( i1d_buffer)))
            time_msec_day = ( mod ( i1d_buffer , microsec_per_day ) ) / 1000

            ! make data missing of missing scan time
            where (i1d_buffer < 0)
               time_msec_day = missing_value_int
            endwhere
            if (.not. allocated ( out % geo % scan_time ) ) allocate (  out % geo % scan_time (dim_seg(2)) )

            out % geo % scan_time =  (/( time_msec_day( ( k -1 ) /16  + 1) , k =ny_start , ny_end )/) 
            deallocate ( time_msec_day )
         case(4) 
            if (.not. allocated ( out % geo % sataz) ) allocate ( out % geo % sataz (dim_seg(1), dim_seg(2)) )
            out % geo % sataz = r2d_buffer 
         case(5)
            if (.not. allocated ( out % geo % satzen) ) allocate ( out % geo % satzen (dim_seg(1), dim_seg(2)) )
            out % geo % satzen = r2d_buffer 
         case(6) 
            if (.not. allocated ( out % geo % solaz) ) allocate ( out % geo % solaz (dim_seg(1), dim_seg(2)) )
            out % geo % solaz = r2d_buffer 
         case(7) 
            if (.not. allocated ( out % geo % solzen) ) allocate ( out % geo % solzen (dim_seg(1), dim_seg(2)) )
            out % geo % solzen = r2d_buffer 
         end select
         
         if (i_gmtco /= 3)  deallocate ( r2d_buffer)
         if (i_gmtco == 3) deallocate ( i1d_buffer)
      end do  
      
      
   
      ! - mbands
      orbit_identifier = trim ( config %file_gmtco_base(11:37))//"*"
      is_mband_on = config %  chan_on_rfl_mband
      
      data_scaled_mband = .true.
      out % file_exists % svm_file_exists (:) = .true.
      
      ! - channel 13 is the only without Factors in VIRRS file
      data_scaled_mband(13) = .false.
      
      do i_mband = 1 , 16
         
         out % mband (i_mband) % is_read = .false.
         
         if ( .not. is_mband_on(i_mband) ) cycle
         
         write ( band_nr_file , '(i2.2)'  )  i_mband
         write ( band_nr_var , '(i2)' )  i_mband
      
         if ( i_mband < 10) write ( band_nr_var , '(i1)' )       i_mband
     
         file_arr_dummy => file_search (trim(config %dir_1b), 'SVM'//trim(band_nr_file)//'_*'//trim(orbit_identifier) , n_files  ) 
      
         if ( n_files == 0 ) out % file_exists % svm_file_exists (i_mband) = .false.
         if ( n_files == 0 ) cycle
         
         file_mband = file_arr_dummy(1)
    
         setname_mband = trim('All_Data/VIIRS-M'//trim(band_nr_var)//'-SDR_All/Reflectance' )
         if ( i_mband > 11) setname_mband = trim('All_Data/VIIRS-M'//trim(band_nr_var)//'-SDR_All/Radiance' )
         setname_mband_fac = trim(setname_mband)//'Factors'
         
         if ( data_scaled_mband ( i_mband ) ) then
            
            call h5readdataset ( trim(config %dir_1b)//trim(file_mband), setname_mband , offset_mband , dim_seg , i2d_buffer  )
            dim_buf = shape ( i2d_buffer )
            allocate ( invalid_pixel ( dim_buf(1) , dim_buf(2) ) )
           
            invalid_pixel = i2d_buffer >= scaled_missing 
            call h5readdataset ( trim(config %dir_1b)//trim(file_mband), setname_mband_fac , factors )
            allocate  ( r2d_data  ( dim_buf(1) , dim_buf(2) )) 
            r2d_data =  i2d_buffer * factors(1) + factors(2)
            deallocate ( i2d_buffer )
         else
            call h5readdataset ( trim(config %dir_1b)//trim(file_mband), setname_mband , offset_mband , dim_seg , r2d_buffer  )
            dim_buf = shape ( r2d_buffer )
            allocate ( invalid_pixel ( dim_buf(1) , dim_buf(2) ) )
            invalid_pixel = r2d_buffer <= unscaled_missing 
            allocate  ( r2d_data  ( dim_buf(1) , dim_buf(2) )) 
            r2d_data = r2d_buffer
            deallocate ( r2d_buffer)
         end if
        
         
         if ( i_mband  <= 11) then 
            if (.not. allocated ( out % mband (i_mband) % ref) ) allocate ( out % mband (i_mband) % ref(dim_seg(1), dim_seg(2)) )
            out % mband ( i_mband) % ref =  100. * r2d_data 
            
            where ( invalid_pixel )  out % mband ( i_mband) % ref  = missing_value_real4
            
         else
            call convert_viirs_radiance ( r2d_data , config % Nu_List(i_mband) , missing_value_real4 )
            if (.not. allocated ( out % mband (i_mband) % rad) ) allocate ( out % mband (i_mband) % rad(dim_seg(1), dim_seg(2)) )
            out % mband ( i_mband) % rad =  r2d_data 
            where ( invalid_pixel ) out % mband ( i_mband) % rad  = missing_value_real4
         end if
         
         if ( allocated ( invalid_pixel) ) deallocate ( invalid_pixel )
         if ( allocated (r2d_data) ) deallocate ( r2d_data )
         
         out % mband (i_mband) % is_read = .true.
      end do
   
      !- ibands
      is_iband_on =config %  chan_on_iband
      dim_seg_iband = dim_seg * 2
      nx_start_iband = 1
      ny_start_iband = ny_start * 2 - 1
      offset_iband = [ nx_start_iband -1 , ny_start_iband - 1 ]
      out % file_exists % svi_file_exists (:) = .true.
      do i_iband = 1 , 5
         if ( .not. is_iband_on(i_iband) ) cycle
         write ( band_nr_file , '(i2.2)'  )  i_iband
         write ( band_nr_var , '(i1)' )       i_iband
         file_arr_dummy => file_search (trim(config %dir_1b), 'SVI'//trim(band_nr_file)//'_*'//trim(orbit_identifier) , n_files  )
         if ( n_files == 0 ) out % file_exists % svi_file_exists (i_iband) = .false.
         if ( n_files == 0 ) cycle
         file_iband = file_arr_dummy(1)

         setname_iband = trim('All_Data/VIIRS-I'//trim(band_nr_var)//'-SDR_All/Reflectance' )
         if ( i_iband > 3) setname_iband = trim('All_Data/VIIRS-I'//trim(band_nr_var)//'-SDR_All/BrightnessTemperature' )
          setname_iband_fac = trim(setname_iband)//'Factors'

         call h5readdataset ( trim(config %dir_1b)//file_iband, setname_iband , offset_iband  , dim_seg_iband , i2d_buffer )
         dim_buf = shape ( i2d_buffer )
         allocate ( invalid_pixel ( dim_buf(1) , dim_buf(2) ) )
         invalid_pixel = i2d_buffer >= scaled_missing
         call h5readdataset ( trim(config %dir_1b)//trim(file_iband), setname_iband_fac , factors )
         allocate  ( r2d_data  ( dim_buf(1) , dim_buf(2) ))
         r2d_data =  i2d_buffer * factors(1) + factors(2)

         if ( i_iband  <= 3) then
            if (.not. allocated ( out % iband (i_iband) % ref) ) allocate ( out % iband (i_iband) % ref(dim_seg_iband(1), dim_seg_iband(2)) )
            out % iband ( i_iband) % ref =  100. * r2d_data
            where ( invalid_pixel )  out % iband ( i_iband) % ref  = missing_value_real4
         else
            if (.not. allocated ( out % iband (i_iband) % bt) ) allocate ( out % iband (i_iband) % bt (dim_seg_iband(1), dim_seg_iband(2)) )
            out % iband ( i_iband) % bt =  r2d_data
            where ( invalid_pixel ) out % iband ( i_iband) % bt  = missing_value_real4
         end if

          deallocate ( invalid_pixel )
          deallocate ( i2d_buffer )
          deallocate ( r2d_data )

          out % iband (i_iband) % is_read = .true.

      end do

    
      ! - dnb
      ! - 
      if (config %  chan_on_dnb ) then
         
         ! - mapping file ( maps from dnb to M-bands resolution)
         file_dnb_idx = trim(config % Ancil_Data_Dir)//'static/viirs/dnb2m_indx.txt'
         lun = getlun()
         dim_seg_dnb(1) = 4064
         dim_seg_dnb(2) = dim_seg(2)
         
         open (unit = lun , file=trim ( file_dnb_idx) , status="old",action="read")
         read (unit = lun , fmt=* , iostat = io_err_stat) d2m_indx
         close (unit = lun)
         
         !- gdnbo products
         out % file_exists % gdnbo_file_exists = .true.
         file_arr_dummy => file_search (trim(config %dir_1b), 'GDNBO*'//trim(orbit_identifier) , n_files )
         if ( n_files > 0 ) then
            file_gdnbo = file_arr_dummy(1)
            call h5readdataset ( trim(config %dir_1b)//file_gdnbo,'All_Data/VIIRS-DNB-GEO_All/LunarAzimuthAngle' ,offset_mband , dim_seg_dnb , r2d_buffer )
            where ( r2d_buffer < unscaled_missing ) r2d_buffer = missing_value_real4
            do i_map = 1, dim_seg(1)
               if (.not. allocated ( out % geo % lunaz) ) allocate ( out % geo % lunaz (dim_seg(1), dim_seg(2)) )
               out % geo % lunaz ( i_map , 1 : ( ny_end - ny_start + 1 ) ) = r2d_buffer ( d2m_indx (i_map) , : )
            end do
            deallocate ( r2d_buffer )
         
            call h5readdataset ( trim(config %dir_1b)//file_gdnbo,'All_Data/VIIRS-DNB-GEO_All/LunarZenithAngle' ,offset_mband , dim_seg_dnb, r2d_buffer )
            where ( r2d_buffer < unscaled_missing ) r2d_buffer = missing_value_real4
            do i_map = 1, dim_seg(1)
               if (.not. allocated ( out % geo % lunzen ) ) allocate ( out % geo %  lunzen (dim_seg(1), dim_seg(2)) )
               out % geo %  lunzen ( i_map , 1 : ( ny_end - ny_start + 1 ) ) = r2d_buffer (d2m_indx ( i_map ) , : )
            end do
           
            deallocate ( r2d_buffer )
        
            allocate ( out % geo % lunrelaz ( dim_seg(1) , dim_seg(2) ) )
            
            call compute_relative_azimuth_viirs ( out % geo % lunaz , out % geo % sataz, out % geo % lunrelaz ) 
           
            call H5ReadDataset( trim(config %dir_1b)//file_gdnbo , 'All_Data/VIIRS-DNB-GEO_All/MoonIllumFraction', out % geo % Moon_Illum_Frac )
            
            call H5ReadDataset( trim(config %dir_1b)//file_gdnbo , 'All_Data/VIIRS-DNB-GEO_All/MoonPhaseAngle', out % geo % Moon_Phase_Angle )
         else
            out % file_exists % gdnbo_file_exists = .false.
         end if

         out % file_exists % svdnb_file_exists = .true.
         file_arr_dummy => file_search (trim(config %dir_1b), 'SVDNB*'//trim(orbit_identifier) , n_files  )
         if ( n_files > 0 ) then
            file_svdnb = file_arr_dummy(1)
        
            call h5readdataset ( trim(config %dir_1b)//file_svdnb, 'All_Data/VIIRS-DNB-SDR_All/Radiance' ,offset_mband , dim_seg_dnb, r2d_buffer )
            where ( r2d_buffer < unscaled_missing ) r2d_buffer = missing_value_real4
            ! --- remap dnb using indexes
            if (.not. allocated ( out % dnb_mgrid % rad ) ) allocate ( out % dnb_mgrid %  rad (dim_seg(1), dim_seg(2)) )
            do i_map = 1, dim_seg(1)
              out % dnb_mgrid % rad(i_map , 1 : ( ny_end - ny_start + 1 ) ) = r2d_buffer(d2m_indx( i_map ) ,  :  )
            end do
            deallocate ( r2d_buffer )
         
            if (.not. allocated ( out % dnb_mgrid % ref ) ) allocate ( out % dnb_mgrid %  ref (dim_seg(1), dim_seg(2)) )
         
            !call read_viirs_date_time(file_gdnbo  , doy = day_of_year  )
            call convert_rad_2_sol_ref_dnb ( out % dnb_mgrid % rad &
                                , out % geo % solzen &
                                , day_of_year &
                                , missing_value_real4 &
                                , out % dnb_mgrid % ref )
         else
            out % file_exists % svdnb_file_exists = .false.
         end if
      end if
      
      ! - cloud products - should be also on a different location!   
      ! - cloud_mask_1b is bad defined, which values can it have sym%YES and NO or 0,1,2?
      ! - goal is to populate cld_mask_aux , cld_phase_aux and cld_type_aux
     
      if ( config %  viirs_cloud_mask_on ) then 
         out % file_exists % iicmo_file_exists = .true.
         file_arr_dummy => file_search (trim(config %dir_1b), 'IICMO*'//trim(orbit_identifier) , n_files  )
         allocate ( out % prd % cld_phase ( dim_seg(1) , dim_seg(2) ) )
         allocate ( out % prd % cld_mask ( dim_seg (1), dim_seg(2) ) )
         allocate ( out % prd % cld_type ( dim_seg(1), dim_seg(2) ) )
         out % prd % cld_mask = missing_value_int
         out % prd % cld_phase = missing_value_int
         out % prd % cld_type = missing_value_int
         
         if ( n_files /= 0 ) then
            file_iicmo = file_arr_dummy(1)
            setname_iicmo='All_Data/VIIRS-CM-IP_All/QF1_VIIRSCMIP'
            call h5readdataset ( trim(config %dir_1b)//file_iicmo, setname_iicmo  , offset_mband , dim_seg , i2d_buffer )
           
            out % prd % cld_mask = 0
 
            ! extract VCM bits
            where ( btest ( i2d_buffer , 2 ) )
               out % prd % cld_mask  =   out % prd % cld_mask  + 1 
            end where
            where  ( btest ( i2d_buffer , 3 ) )
               out % prd % cld_mask  =   out % prd % cld_mask  + 2 
            end where

            ! assign missing values
            where ( i2d_buffer == 0 )
               out % prd % cld_mask = missing_value_int
            end where

            deallocate ( i2d_buffer )
            
            setname_iicmo='All_Data/VIIRS-CM-IP_All/QF6_VIIRSCMIP'
            call h5readdataset ( trim(config %dir_1b)//file_iicmo, setname_iicmo , offset_mband , dim_seg , i2d_buffer )
            shape_buffer = shape (i2d_buffer)
            allocate ( cld_type_idps ( shape_buffer(1), shape_buffer(2)) )
            
            ! - cascade to avoid  seg fault due to too high memory 
            cld_type_idps = log2int ( btest( i2d_buffer , 0 ) )
            cld_type_idps =cld_type_idps +  2.*  log2int(btest( i2d_buffer , 1 ) )
            cld_type_idps =cld_type_idps +  4.*  log2int(btest( i2d_buffer , 2 ) )
         
           ! cld_type_idps = log2int ( btest( i2d_buffer , 0 ) ) &
           !                    + 2.*  log2int(btest( i2d_buffer , 1 ) ) &
           !                    + 4. * log2int(btest( i2d_buffer , 2 ) )
            
            out % prd % cld_phase  = 5
            where ( cld_type_idps == 1 .or. cld_type_idps == 2 )
               out % prd % cld_phase  = 0
            end where
        
            where ( cld_type_idps == 3 )
               out % prd % cld_phase  = 1
            end where
        
            where ( cld_type_idps == 4 )
               out % prd % cld_phase  = 2
            end where
         
            where ( cld_type_idps == 5 .or. cld_type_idps == 6 .or. cld_type_idps == 7 )
               out % prd % cld_phase  = 4
            end where

            ! assign missing values
            where ( i2d_buffer == 0 )
               out % prd % cld_phase = missing_value_int
            end where

            ! - type
      
            out % prd % cld_type  = 10
            
            where ( cld_type_idps == 1 )
               out % prd % cld_type = 0
            end where     
            
            where ( cld_type_idps == 2 )
               out % prd % cld_type = 1
            end where
            
            where ( cld_type_idps == 3 )
               out % prd % cld_type = 3
            end where
            
            where ( cld_type_idps == 4 )
               out % prd % cld_type = 4
            end where
            
            where ( cld_type_idps == 5 )
               out % prd % cld_type = 6
            end where
            
            where ( cld_type_idps == 6 )
               out % prd % cld_type = 7
            end where
            
            where ( cld_type_idps == 7 )
               out % prd % cld_type = 8
            end where
            
            deallocate ( cld_type_idps )
 
            where  ( out % prd % cld_type == 10 )
               out % prd % cld_mask = missing_value_int
            end where
            
            ! assign missing values
            where ( i2d_buffer == 0 )
               out % prd % cld_type = missing_value_int
            end where

            deallocate ( i2d_buffer )
            
         else 
            out % file_exists % iicmo_file_exists = .false.
            print*,'IICMO file not found , cld_mask_aux and cld_type_aux are set to missing'
            
         end if  
      end if
      ! -horrible global parameters this has to be gone soon! 
      !Num_Scans_Read = Ny_End - Ny_Start + 1
      ! Fill in Scan_Number
      !do i = 1, num_scans_per_segment
      !   scan_number(i) = ny_start + i - 1
      !end do
     
      if ( .not. gap_pattern_computed  ) then
         call compute_viirs_bowtie_gap_pattern()
         gap_pattern_computed = .true.
      end if   
   
      call fill_viirs_bowtie_gaps  ( ny_start ,  dim_seg(2) , out )
     
      
      ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !  everything between these xxx have to go to other locations  
            

      
      !- ascending  (global varaibel )
      allocate ( out % geo % ascend ( dim_seg (2) ) ) 
      out % geo % ascend = 0  
      do i = 1 , ny_end -ny_start - 1
         if ( out%geo%lat (dim_seg(1) / 2 , i + 1) <= out%geo%lat(dim_seg(1) / 2 , i ) ) out % geo % ascend  ( i )  = 1
      end do

      allocate ( out % geo % relaz ( dim_seg(1) , dim_seg(2) ) )
    
      ! rel azimuths  - these are all global variables
      call  compute_relative_azimuth_viirs ( out%geo%solaz , out%geo%sataz , out%geo%relaz)
     
   end subroutine read_viirs_level1b
   
   
   !-------------------------------------------------------------------------------------
   ! subroutine to compute the bowtie gap pattern.  
   ! This pattern repeats every 48 scans 
   ! all VIIRS files should be integer multiples of this pattern
   !
   ! Gap_Pixel_Mask_Pattern = a binary mask that identifies these bowtie gaps
   ! Gap_Line_Idx = line index for each pixel in pattern including gap pixels
   !
   ! A description of the pattern
   !
   !--------  line type 3
   !----      line type 4
   ! 12 lines without gaps
   !----      line type 1
   !--------  line type 2
   !--------  line type 3
   !----      line type 4
   ! 12 lines without gaps
   !----      line type 1
   !--------  line type 2
   !--------  line type 3
   !----      line type 4
   ! 12 lines without gaps
   !----      line type 1
   !--------  line type 2
   !
   !  line types 1 and 4 have 1280 missing pixels
   !  line types 2 and 3 have 2016 missing pixels
   ! 
   !-------------------------------------------------------------------------------------
   subroutine COMPUTE_VIIRS_BOWTIE_GAP_PATTERN()
    
      integer (kind=int4), dimension(Ny_Pattern):: Line_Type
      integer (kind=int4), parameter:: Ngap_1 = 640  !1280
      integer (kind=int4), parameter:: Ngap_2 = 1008 !2016
      integer (kind=int4), parameter:: Ngap_3 = 1008 !2016
      integer (kind=int4), parameter:: Ngap_4 = 640  !1280
      integer (kind=int4):: iline
      integer (kind=int4):: i1
      integer (kind=int4):: i2

      !--- define the line patterns as described above
      Line_Type = (/3,4,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,1,2/)

  
      if (.not. allocated(Gap_Line_Idx_Pattern)) Allocate(Gap_Line_Idx_Pattern(Nx_Pattern,Ny_Pattern))    
      if (.not. allocated(Gap_Pixel_Mask_Pattern)) Allocate(Gap_Pixel_Mask_Pattern(Nx_Pattern,Ny_Pattern))    

      do iline = 1 , Ny_Pattern

         Gap_Line_Idx_Pattern( : , iline ) =  -999 
         Gap_Pixel_Mask_Pattern( : , iline ) = .false.

         if (line_Type(iline) == 1) then
            
            i1 = 1
            i2 = Ngap_1 
            Gap_Line_Idx_Pattern( i1 : i2 , iline) = iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , iline) = .true.

            i1 = Nx_Pattern - Ngap_1 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , iline ) = iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , iline ) = .true.
         end if

         if (Line_Type(iline) == 2) then 
            i1 = 1
            i2 = Ngap_1
            Gap_Line_Idx_Pattern( i1 : i2 , iline ) = iline - 2
            Gap_Pixel_Mask_Pattern( i1 : i2 , iline ) = .true.

            i1 = Ngap_1 + 1
            i2 = Ngap_2
            Gap_Line_Idx_Pattern( i1 : i2 , iline ) = iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , iline ) = .true.

            i1 = Nx_Pattern - Ngap_1 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , iline ) = iline - 2
            Gap_Pixel_Mask_Pattern( i1 : i2 , iline ) = .true.

            i1 = Nx_Pattern - Ngap_2 + 1
            i2 = i1  + (Ngap_2 - Ngap_1)
            Gap_Line_Idx_Pattern( i1 : i2 , iline ) = iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , iline ) = .true.
         end if

         if (Line_Type(iline) == 3) then
            i1 = 1
            i2 = Ngap_4
            Gap_Line_Idx_Pattern(i1:i2,iline) = iline + 2
            Gap_Pixel_Mask_Pattern(i1:i2,iline) = .true.

            i1 = Ngap_4 + 1
            i2 = Ngap_3
            Gap_Line_Idx_Pattern(i1:i2,iline) = iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,iline) = .true.

            i1 = Nx_Pattern - Ngap_4 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern(i1:i2,iline) = iline + 2
            Gap_Pixel_Mask_Pattern(i1:i2,iline) = .true.

            i1 = Nx_Pattern - Ngap_3 + 1
            i2 = i1 + (Ngap_3 - Ngap_4)
            Gap_Line_Idx_Pattern(i1:i2,iline) = iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,iline) = .true.
         end if

         if (Line_Type(iline) ==  4) then
            i1 = 1
            i2 = Ngap_4
            Gap_Line_Idx_Pattern(i1:i2,iline) = iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,iline) = .true.

            i1 = Nx_Pattern - Ngap_4 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , iline ) = iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,iline) = .true.
         end if
       
      end do    

   end subroutine compute_viirs_bowtie_gap_pattern
!------------------------------------------------------------------------------
! this routine uses the bowtie gap pattern, applies it to an arbitrary 
! segment of data and fills in the observations with the closest valid data
!------------------------------------------------------------------------------
    subroutine fill_viirs_bowtie_gaps ( Line_Start, Number_of_Lines  , out )
      integer(kind=int4), intent(in):: Line_Start
      integer(kind=int4), intent(in):: Number_of_Lines
      type (  viirs_data_out ) , intent (inout) :: out
   
      integer(kind=int4):: Line_Offset
      integer(kind=int4):: Line_in_Pattern
      integer(kind=int4):: Line_in_Segment
      integer(kind=int4):: Line_Idx
      integer(kind=int4):: Elem_Idx
      integer :: Num_Pix = 3200

      integer , dimension(:,: ) , allocatable  :: gap_line_idx
      logical , dimension(:,: ) , allocatable  :: gap_pixel_mask
     
      integer :: i_mband
      integer :: missing_value_int1 = -999
    
      allocate ( gap_pixel_mask (3200, number_of_lines) , gap_line_idx(3200, number_of_lines))
      
      do line_in_segment = 1 ,  number_of_lines
         
         line_in_pattern = mod (line_start -1 + line_in_segment , ny_pattern)
         if (Line_In_Pattern == 0) Line_In_Pattern = Ny_Pattern
      
         line_offset = line_in_segment - line_in_pattern
        
         gap_line_idx ( : , line_in_segment ) = Gap_Line_Idx_Pattern(:,Line_in_Pattern) + line_offset
         
         where ( gap_line_idx(:,Line_in_Segment) <= 0 )
            Gap_Line_Idx(:,Line_in_Segment) = 1
         end where
         
         where ( gap_line_idx(:,Line_in_Segment) > Number_of_Lines )
            Gap_Line_Idx(:,Line_in_Segment) = Number_of_Lines
         end where
         
         ! - write gap mask to output
         Gap_Pixel_Mask(:,Line_in_Segment) = Gap_Pixel_Mask_Pattern(:,Line_in_Pattern)
      
         do Elem_Idx = 1, Num_Pix
            if (Gap_Pixel_Mask(Elem_Idx,Line_in_Segment) ) then
               Line_Idx = Gap_Line_Idx(Elem_Idx,Line_in_Segment)
               
               do i_mband = 1 , 11
                  if ( .not. allocated( out % mband(i_mband) % ref ) &
                       .or.  .not. out % file_exists % svm_file_exists (i_mband)   ) cycle
                  out % mband (i_mband) % ref (elem_idx , line_in_segment) = out%mband(i_mband)%ref(elem_idx,line_idx)
               end do 
               
               do i_mband = 12,16 
                  if ( .not. allocated (out % mband(i_mband) % rad ) &
                       .or. .not. out % file_exists % svm_file_exists (i_mband) ) cycle
                  out % mband(i_mband) % rad (elem_idx,line_in_segment) = out % mband(i_mband)%rad(elem_idx,line_idx)
               end do
               
               if ( allocated (out % prd % cld_mask) .and. out % file_exists % iicmo_file_exists ) then 
                  out % prd % cld_type(elem_idx,line_in_segment) = out % prd % cld_type(elem_idx,line_idx)
                  out % prd % cld_phase(elem_idx,line_in_segment) = out % prd % cld_phase(elem_idx,line_idx)
                  out % prd % cld_mask(Elem_Idx,Line_in_Segment) = out % prd %  cld_mask(Elem_Idx,Line_Idx)
               end if
               
            end if
         end do
      
      end do
      
      
      if ( allocated ( out % prd % cld_mask ) .and. out % file_exists % iicmo_file_exists  ) then
         where( out % prd % cld_mask == Missing_Value_Int1) 
            Gap_Pixel_Mask = .true.
         end where
      end if
      
      deallocate (  gap_line_idx)
       allocate ( out % gap % mask (3200, number_of_lines))
      out % gap % mask = gap_pixel_mask
      deallocate (  gap_pixel_mask)
      
   end subroutine fill_viirs_bowtie_gaps
 
 
   
   !  this routine should be at a different place
   subroutine  compute_relative_azimuth_viirs ( ang1 , ang2, rel_az_out )
      real , dimension(:,:), intent(in) :: ang1
      real , dimension(:,:), intent(in) :: ang2
      real , dimension(:,:)  :: rel_az_out
    
    
      rel_az_out = abs (ang1 - ang2 )
      where ( rel_az_out > 180 ) rel_az_out = 180.0 - rel_az_out
     
   end subroutine compute_relative_azimuth_viirs
  
 !-----------------------------------------------------------------------------------------
   !  Extract time information from VIIRS filename - should explore use of header for this
   !   assumingly called from outside  ...
   !-----------------------------------------------------------------------------------------
  
   subroutine READ_VIIRS_DATE_TIME_ATT (Path, Infile, Year , Doy , Start_Time &
                , End_Time , Orbit , Orbit_Identifier , End_Year, End_Doy )
      ! Get the date & time from the file's name
      implicit none
      
      character(len=*), intent(in) :: Path
      character(len=*), intent(in) :: Infile   
      integer, intent(out) , optional :: Year
      integer, intent(out)  , optional:: Doy    !day of year
      integer, intent(out) , optional :: Start_Time  !millisec
      integer, intent(out)  , optional:: End_Time    !millisec
      integer, intent(out)  , optional:: Orbit
      character(38), intent(out) , optional :: Orbit_Identifier
      integer , intent(out) , optional :: End_Year
      integer, intent(out)  , optional:: End_Doy    !day of year
  
      integer :: Month
      integer :: Day
      integer :: Start_Hour
      integer :: Start_Minute
      integer :: Start_Sec

      integer :: Days_Of_Year
      integer :: End_Hour
      integer :: End_Minute
      integer :: End_Sec
      integer :: Year_Loc
      integer :: Doy_Loc    !day of year
      integer :: End_Year_Loc
      integer :: End_Doy_Loc    !day of year
      integer :: Start_Time_Loc  !millisec
      integer :: End_Time_Loc    !millisec
      integer :: Orbit_Loc
      character(38):: Orbit_Identifier_Loc
      character(50) :: String_Tmp
 

      print *,trim(Path)
      print *,trim(Infile)

      ! --- read time and date from the attributes      
      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Aggr/AggregateBeginningOrbitNumber', Orbit_Loc)

      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Aggr/AggregateBeginningTime', String_Tmp)

      read(String_Tmp(1:2), fmt="(I2)") Start_Hour
      read(String_Tmp(3:4), fmt="(I2)") Start_Minute
      read(String_Tmp(5:6), fmt="(I2)") Start_Sec

      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Aggr/AggregateEndingTime', String_Tmp)

      read(String_Tmp(1:2), fmt="(I2)") End_Hour
      read(String_Tmp(3:4), fmt="(I2)") End_Minute
      read(String_Tmp(5:6), fmt="(I2)") End_Sec

      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/Beginning_Date', String_Tmp)

      read(String_Tmp(1:4), fmt="(I4)") Year_Loc
      read(String_Tmp(5:6), fmt="(I2)") Month
      read(String_Tmp(7:8), fmt="(I2)") Day
      
  
      !         1         2	    3	      4	        5	  6	    7	      8
      !12345678901234567890123456789012345678901234567890123456789012345678901234567890
      !GMODO_npp_d20100906_t2110510_e2112156_b00012_c20110707160532497848_noaa_ops.h5

      !---- store orbit number
      Orbit_Identifier_Loc = Infile(7:44)

      !--- compute day of year
      call JULIAN ( Day, Month, Year_Loc, Doy_Loc )
 

      ! --- Calculate start and end time
      Start_Time_Loc = ((Start_Hour * 60 + Start_Minute) * 60 + Start_Sec) * 1000
      End_Time_Loc = ((End_Hour * 60 + End_Minute) * 60 + End_Sec) * 1000
  
      End_Doy_Loc = Doy_Loc
      End_Year_Loc = Year_Loc
      if ( End_time_loc <= start_time_loc) then
         End_Doy_Loc = End_Doy_Loc + 1
         Days_Of_Year = 365
         if ( modulo(Year_Loc,4) == 0)  Days_Of_Year = 366
         if ( End_Doy_Loc > Days_Of_Year) then
            End_Doy_Loc = 1
            End_Year_Loc = End_Year_Loc + 1
         end if  
      end if
      
      if ( present ( Year )) Year = Year_Loc
      if ( present ( Doy )) Doy = Doy_Loc
      if ( present ( Start_Time )) Start_Time = Start_Time_Loc
      if ( present ( End_Time )) End_Time = End_Time_Loc
      if ( present ( Orbit )) Orbit = Orbit_Loc
      if ( present ( Orbit_Identifier )) Orbit_Identifier = Orbit_Identifier_Loc
      if ( present ( End_Year )) End_Year = End_Year_Loc
      if ( present ( End_Doy )) End_Doy = End_Doy_Loc
      
   end subroutine READ_VIIRS_DATE_TIME_ATT


!---------------------------------------------------------------------------------
!  subroutine READ_NUMBER_OF_SCANS_FROM_VIIRS ( Infile, Number_Of_Viirs_Lines, Error_Out )
!  to read number of scans,  called from the bridge
!---------------------------------------------------------------------------------
   SUBROUTINE READ_NUMBER_OF_SCANS_FROM_VIIRS ( Infile, Number_Of_Viirs_Lines, Error_Out )
   
      CHARACTER(Len=*), INTENT(IN) :: Infile  
      INTEGER(kind=int4), INTENT(OUT) :: Error_Out
      INTEGER(KIND=INT4), INTENT(OUT):: Number_of_Viirs_Lines
      CHARACTER(Len=100) :: Setname
      integer ,dimension(:), pointer ::test


      error_out = 0
      Setname = 'All_Data/VIIRS-MOD-GEO-TC_All/NumberOfScans'
      call H5ReadDataset( infile, setname, test )
     
      Number_of_Viirs_Lines = sum(test)
   
   END SUBROUTINE READ_NUMBER_OF_SCANS_FROM_VIIRS

   !====================================================================
   ! Function Name: CONVERT_VIIRS_RADIANCE
   !
   ! Function:
   !    Convert to units of the VIIRS radiance values from the that used
   !    in the IDPS level-1b to that expected by CLAVR-x
   !
   ! Description: 
   !   
   ! Calling Sequence: rad_new = convert_viirs_radiance(rad_old,nu,missing_value)
   !   
   !
   ! Inputs:
   !   rad_old = radiance in units of W/m^2/micron/str (2d array)
   !   nu = channels equivalent width in units of cm^-1
   !   missing_value = value assigned to missing radiance values
   !
   ! Outputs: 
   !   rad_new = radiance in units of mW/m^2/cm^-1/str (2d array)
   !
   ! Dependencies:
   !
   ! Restrictions:  None
   !
   ! Reference: algebraic manipulation of Planck Equation
   ! ---------------------------------------------------------------------------------------
   subroutine CONVERT_VIIRS_RADIANCE(radiance,nu,missing_value)
      real (kind=real4), dimension(:,:), intent(inout):: Radiance
      real (kind=real4), intent(in):: Nu
      real (kind=real4), intent(in):: Missing_Value

      where(Radiance /= Missing_Value)
         Radiance = Radiance * (((10000.0 / Nu )**2) / 10.0)
      end where

      return

   end subroutine CONVERT_VIIRS_RADIANCE



    subroutine convert_rad_2_sol_ref_dnb ( radiance , solzen , day_of_year , missing_value , reflectance )
   
      real (kind=real4), dimension(:,:), intent(in) :: radiance
      real (kind=real4), dimension(:,:), intent(in) :: solzen
      integer , intent(in) :: day_of_year
      real (kind=real4), intent(in) :: missing_value
      real (kind=real4), dimension(:,:), intent(out) :: reflectance
      real (kind=real4), parameter :: fo = 0.044217282   ! dnb solar energy in w/cm2 ( Source? )
      real , parameter :: PI = 3.14159265359
      real , parameter :: DTOR = PI / 180.
      real :: sun_earth_distance
     
      sun_earth_distance = 1.0 - 0.016729 * cos ( 0.9856 * ( day_of_year - 4.0 ) * DTOR )

      reflectance = missing_value

      where(radiance /= missing_value)
!         reflectance = 100*(pi*radiance*sun_earth_distance**2) / (cosd(solzen)*fo)
         reflectance = 100*(pi*radiance*sun_earth_distance**2) / (cos(solzen*pi/180.0)*fo)
      end where

   end subroutine convert_rad_2_sol_ref_dnb

  !---------------------------------------------------------------------------------------
  !    transforms LOGICAL data type to an integer
  !
   function log2int (b) result (r)
      logical , dimension(:,:) :: b
      integer , dimension(:,:),allocatable :: r
      integer :: dims(2)
    
      dims = shape ( b)
      allocate ( r (dims(1), dims(2) ) )
      r = 0
      where ( b )
         r = 1
      end where 
   end function log2int

!---------------------------------------------------------------------------------------

!-------------------------------------------------
! subroutine JULIAN(iday,imonth,iyear,jday)
! compute julian day
! input:
!         iday - integer day
!         imonth - integer month
!         iyear - integer year (2 or four digits)
!         
! output : jday - julian day
!--------------------------------------------------
 subroutine JULIAN(iday,imonth,iyear,jday)

!-- Computes julian day (1-365/366)
        integer, intent(in)::  iday,imonth,iyear
        integer, intent(out):: jday
        integer::  j
        integer, dimension(12)::  jmonth

        jmonth = reshape ((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))

        jday = iday
        if (modulo(iyear,4) == 0) then
            jmonth(2)=29
        endif

        do j = 1,imonth-1
           jday = jday + jmonth(j)
        end do

   end subroutine JULIAN
   !
   !  public routine to deallocate output structure
   !  history: 04/30/2013 AW
   ! 
   subroutine dealloc_viirs_data_out ( this )
      class ( viirs_data_out) :: this
      integer :: m 
      
      
      if ( allocated (this%prd%cld_mask) ) deallocate (this%prd%cld_mask )
      if ( allocated (this%prd%cld_phase) ) deallocate (this%prd%cld_phase )
      if ( allocated (this%prd%cld_type) ) deallocate (this%prd%cld_type )
      
      if (allocated ( this%geo%solzen) ) deallocate ( this%geo%solzen) 
      if (allocated ( this%geo%satzen) ) deallocate ( this%geo%satzen) 
      if (allocated ( this%geo%solaz) ) deallocate ( this%geo%solaz) 
      if (allocated ( this%geo%sataz) ) deallocate ( this%geo%sataz) 
      
      if (allocated ( this%geo%relaz) ) deallocate ( this%geo%relaz) 
      if (allocated ( this%geo%lunaz) ) deallocate ( this%geo%lunaz) 
      if (allocated ( this%geo%lunzen) ) deallocate ( this%geo%lunzen) 
      if (allocated ( this%geo%lunrelaz) ) deallocate ( this%geo%lunrelaz) 
      if (allocated ( this%geo%lat) ) deallocate ( this%geo%lat) 
      if (allocated ( this%geo%lon) ) deallocate ( this%geo%lon) 
      
      if (allocated ( this%geo%scan_time) ) deallocate ( this%geo%scan_time) 
      if (allocated ( this%geo%ascend) ) deallocate ( this%geo%ascend)
      
      if (allocated ( this%dnb%ref_lun)) deallocate ( this%dnb%ref_lun)
      if (allocated ( this%dnb%ref_sol)) deallocate ( this%dnb%ref_sol)
      if (allocated ( this%dnb%rad)) deallocate ( this%dnb%rad)
      
       
      if (allocated ( this%dnb_mgrid%ref)) deallocate ( this%dnb_mgrid%ref)
      if (allocated ( this%dnb_mgrid%rad)) deallocate ( this%dnb_mgrid%rad)
      
      if ( allocated (this%gap%mask) ) deallocate (this%gap%mask )
      
      do m = 1 , 16
        if (allocated (this%mband (m) %ref )) deallocate ( this%mband (m) %ref )
        if (allocated (this%mband (m) %rad )) deallocate ( this%mband (m) %rad )
        if (allocated (this%mband (m) %bt )) deallocate ( this%mband (m) %bt )
      end do  
      
       do m = 1 , 5
        if (allocated (this%iband (m) %ref )) deallocate ( this%iband (m) %ref )        
        if (allocated (this%iband (m) %bt )) deallocate ( this%iband (m) %bt )
      end do  
     
   end subroutine dealloc_viirs_data_out
   
   
   
   subroutine dealloc_mband_str ( this )
      class ( mband_str ) :: this
      print*,'mband deallocation'
      if ( allocated (this%ref) ) deallocate (this%ref )
      if ( allocated (this%rad) ) deallocate (this%rad )
      if ( allocated (this%bt) ) deallocate (this%bt )
      print*,'done'
   end subroutine dealloc_mband_str 

end module viirs_read_mod

