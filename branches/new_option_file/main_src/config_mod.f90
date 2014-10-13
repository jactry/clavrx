! $Id:$
module config_mod
   
  
  ! - AVHRR specific 
   type conv_avhhr_type
      logical :: active_rfl_cal_1b
      logical :: use_therm_cal_1b
      logical :: fill_bytes_1b
   end type conv_avhhr_type
   
   type conf_process_type
      logical :: run_sst
      logical :: run_cld
      logical :: run_aot
      logical :: run_erb
      logical :: run_ash
   end type conf_process_type
   
   type conf_diag_type
      logical :: is_diag_on
      logical :: is_node_des
      real :: min_lat
      real :: max_lat 
   end type conf_diag_type
   
   type conf_sfc_type
      integer :: oisst_mode
      logical :: use_seebor
      logical :: use_hres_cloudtype
      logical :: use_hres_land_mask
      logical :: use_hres_coast_mask
      logical :: use_hres_elev
      logical :: use_volc_mask 
      logical :: use_dark_compo 
      integer :: snow_mask_mode
   end type conf_sfc_type
   
   type conf_mask_type
      logical :: use_modis_clear
      logical :: use_prob_clear_res
      logical :: use_lrc
      logical :: read_aux
      character ( len =256) :: bayesian_mask_classifier
   end type conf_mask_type
   
    type conf_channel_type
      logical:: is_on 
   end type conf_channel_type
   
   type conf_file_type
      character (len=256) :: l1b_path
      character (len=256) :: out_path
      character (len=256) , allocatable, dimension (:)  :: infile
      integer , allocatable, dimension(:) :: wmo_id
      integer, allocatable, dimension (:)  :: ETsensor 
      integer, allocatable, dimension (:)  :: ETresolu 
      character(len=20), allocatable :: sensor_name(:)
   end  type conf_file_type
   
   type conf_algo_modes_type
      integer :: dcomp_mode
      integer :: mask_mode
      integer :: acha_mode
      integer :: nlcomp_mode
   end type conf_algo_modes_type
     
   type conf_user_opt_type
      integer :: n_files
      character (len=256) :: ancil_path
      character (len=256) :: temp_path
      type (conf_algo_modes_type ) :: user_modes
      type (conf_algo_modes_type ) :: updated_modes
      integer :: rtm_mode
      integer :: nwp_mode
      logical :: do_smooth_nwp
      logical :: do_process_cloudfree
      logical :: do_compress_out
      type (conf_process_type ) :: algo_group
      type (conf_diag_type ) :: diag
      type ( conf_sfc_type ) :: sfc
      type ( conf_mask_type) :: mask
      type ( conv_avhhr_type ) :: avhrr
      type ( conf_channel_type ) , dimension(40) :: chan
      logical :: use_gzip     
      logical :: subset_pixel_hdf
      
      type ( conf_file_type ):: file
      
      contains
      procedure, public :: set_config
      procedure, private :: read_user_configuration
      procedure, private :: read_input_files
      procedure, public :: update_modes
      procedure, public :: print_it
   end type conf_user_opt_type
   
   enum , bind(C)
      enumerator :: ETsensor_invalid = 0
      enumerator :: ETsensor_avhrr  = 1
      enumerator :: ETsensor_aqua 
      enumerator :: ETsensor_terra
      enumerator :: ETsensor_seviri 
      enumerator :: ETsensor_goes_first = 5
      enumerator :: ETsensor_goes_mop = 5
      enumerator :: ETsensor_goes_il = 6
      enumerator :: ETsensor_goes_last = 6
      enumerator :: ETsensor_viirs 
      enumerator :: ETsensor_mtsat
      enumerator :: ETsensor_aqua_mac 
      enumerator :: ETsensor_iff_viirs
      enumerator :: ETsensor_iff_aqua
   end enum 
   
   enum, bind(C)
      enumerator :: ETresolu_invalid = 0
      enumerator :: ETresolu_1km
      enumerator :: ETresolu_5km
   
   end enum
   
   
contains
   
   !
   !
   !
   subroutine read_input_files (conf)
      use file_tools
      class ( conf_user_opt_type ) :: conf
      character ( len =256 ) :: input_file = 'file_list'
      integer :: n_files
      integer :: i_file
      integer :: ios
      
     
      n_files = file_nr_lines (trim(input_file)) - 2 
      conf % n_files = n_files
      allocate ( conf % file % infile (n_files) ) 
      allocate ( conf % file % ETsensor ( n_files ) )
      allocate ( conf % file % sensor_name ( n_files ) )
      allocate ( conf % file % wmo_id ( n_files ) )
      
       lun = getlun()
      open(unit=lun,file = trim(input_file ),status="old",action="read",iostat=ios)
      read(unit=lun,fmt="(a)")  conf % file % l1b_path
      read(unit=lun,fmt="(a)")  conf % file % out_path
      
      do i_file =1 , n_files
         read(unit=lun,fmt="(a)")  conf % file % infile (i_file)
         conf % file % ETsensor (i_file) &
                 &  =  sensor_from_filename (trim(conf % file % infile (i_file)) , conf % file % l1b_path )
         conf % file % sensor_name (i_file) = sensorname ( conf % file % ETsensor (i_file) )
      end do
      
      close ( unit = lun )
      
      
   end subroutine read_input_files 
   
   function sensorname (et_sensor) result ( sensorn)
      integer, intent(in) :: et_sensor
      character(len=20) :: sensorn 
      
      select case (et_sensor)
           
      case(ETsensor_invalid)
         sensorn='INVALID'
      case(ETsensor_avhrr)
         sensorn='AVHRR'
      case(ETsensor_aqua)
         sensorn='MODIS-AQUA'
      case(ETsensor_terra)
         sensorn='MODIS-TERRA'
      case(ETsensor_seviri) 
         sensorn='MSG'
      
      case(ETsensor_goes_mop)
         sensorn='GOES-MOP'
      case(ETsensor_goes_il)
         sensorn='GOES-IL'
     
      case(ETsensor_viirs )
         sensorn='VIIRS'
      case(ETsensor_mtsat)
         sensorn='MTSAT'
      case(ETsensor_aqua_mac) 
         sensorn='AQUA_MAC'
      case(ETsensor_iff_viirs)
         sensorn='VIIRS-IFF'
      case(ETsensor_iff_aqua)
         sensorn='AQUA-IFF'
      end select 
   
   end function sensorname 
   
   
   !
   !  determine sensor 
   !
   function sensor_from_filename (file , path ) result (sensor)  
      use area_file_mod
      character (len =* ) , intent(in) :: file
      character (len =* ) , intent(in) :: path
      integer  :: sensor
      integer :: resolu
      integer :: sat_id_num
      
      character ( len =1000) :: file_unzipped
      character ( len =100) :: dir_in
           
      sensor = ETsensor_invalid
      resolu = ETresolu_invalid
      
      !-  test for modis
      if ( index( file , 'MYD') ) sensor = ETsensor_aqua
      if ( index( file , 'MOD') ) sensor = ETsensor_terra
      if ( index( file , 'MAC02') ) sensor = ETsensor_aqua_mac
      
      if ( index( file , 'D02SSH') ) resolu = ETresolu_5km
      if ( index( file , 'D021KM') ) resolu = ETresolu_1km
      
      !- test for VIIRS
      if ( index( file , 'GMTCO') ) sensor = ETsensor_viirs
      
      
      ! - test for PEATE IFF VIIRS/MODIS data set
      if ( index( file , 'IFF_npp') ) sensor = ETsensor_iff_viirs
      if ( index( file , 'IFF_aqua') ) sensor = ETsensor_iff_aqua
      
      if ( sensor /= ETsensor_invalid) return
      
      call uncompress_file ( trim(file) ,   file_unzipped , dir_in = trim(path) ,  tmp_dir_in = 'temp/')
      
      ! - test if this is an area file
      if (   is_area_file ( trim(file_unzipped) , sat_id_num ) ) then
         
         select case (sat_id_num )
         case(51:53)
            sensor = ETsensor_seviri
         case(84,85)
            sensor = ETsensor_mtsat
         case(37,38)
            sensor = ETsensor_fy2
         case(250)
            sensor = ETsensor_coms
         case (71,73,75,77,79,181,183,185)
            sensor = ETsensor_goes_sounder
         case (70,72,74,76)
            sensor = ETsensor_goes_il
         case ( 78, 180,182,184)
            sensor = ETsensor_goes_mop
         end select
      
      end if
      
     
      
     
   end function sensor_from_filename
   
   ! called from process_clavrx
   !
   !
   subroutine set_config ( conf)
      class ( conf_user_opt_type ) :: conf
      
      call conf % read_user_configuration ()
      call conf % read_input_files
     
      
   end subroutine set_config 
   
   !
   !   copy user options to clavrx options
   !  and check if these are possible for this 
   !  sensor and file
   !
   subroutine update_modes ( this, i_file )
      class ( conf_user_opt_type ) :: this
      integer , intent(in) :: i_file
      integer  :: acha , dcomp
      
      this % updated_modes = this % user_modes
      acha = this % user_modes % acha_mode
      dcomp = this % user_modes % dcomp_mode
      
      !------------------------------------------------------------------------
      !--- ACHA MODE Check
      !         (0=none; 1 = 11/12; 2 = 11/13.3; 
      !---       3 = 11,12,13; 4=8.5/11/12; 
      !---       5=11/12/6.7; 6=11/13.3/6.7; 7=11/6.7 8 =11 )
      !  -- dcomp 
      !------------------------------------------------------------------------
      
      select case ( this % file % ETsensor(i_file) )
      case ( ETsensor_avhrr )
       
      case ( ETsensor_goes_il )
          if ( acha == 2 .or. acha == 3 .or. acha == 4 .or. acha == 6 ) then   
            print*, 'acha mode ',acha,' incompatible with GOES-IL'
            print*,'changing to default for GOES-IL'
            this % updated_modes % acha_mode = 5
         end if
         
          if ( dcomp /= 0 .and. dcomp /= 3 ) then
            print*, 'dcomp mode ',dcomp,' incompatible with GOES'
            print*,'changing to default for GOES'
            this % updated_modes % dcomp_mode = 3
         end if
         
      case ( ETsensor_goes_mop )
          if ( acha == 1 .or. acha == 3 .or. acha == 4 .or. acha == 5 ) then   
            print*, 'acha mode ',acha,' incompatible with GOES-MOP'
            print*,'changing to default for GOES-MOP'
            this % updated_modes % acha_mode = 6
         end if
         
         if ( dcomp /= 0 .and. dcomp /= 3 ) then
            print*, 'dcomp mode ',dcomp,' incompatible with GOES'
            print*,'changing to default for GOES'
            this % updated_modes % dcomp_mode = 3
         end if
         
      case ( ETsensor_viirs )
         if ( acha == 2 .or. acha == 3 .or. acha == 5 .or. acha == 6 .or. acha == 7) then   
            print*, 'acha mode ',acha,' incompatible with VIIRS'
            print*,'changing to default for VIIRS'
            this % updated_modes % acha_mode = 4
         end if
      case ( ETsensor_mtsat )
      
      case default
         print*, 'somrhong wrong'
      
      end select 
      
   end subroutine update_modes
   
  
   subroutine read_user_configuration ( conf )
      class ( conf_user_opt_type ) :: conf
      character ( len =256 ) :: option_file = 'clavrx_options'
      integer :: lun
      integer :: ios0
      integer :: errstat
      integer :: int_dummy
      
      ! - check if file is existing
      
      
      ! -- lun - get_lun()
      lun =13
      
      open ( unit = lun , file = trim(option_file) &
                  , iostat = ios0 &
                  , action = 'read' &
                  , status = 'old' &
                  , position = 'rewind' )
                  
      errstat = 0
      
      if ( ios0 /= 0 ) then            
         errstat = 1               
         print*, 'Error opening clavrx_options file'
         stop 1
      else
         read(unit=lun,fmt="(a)") conf % ancil_path
         read(unit=lun,fmt="(a)") conf % temp_path
         read(unit=lun,fmt=*) conf %user_modes% mask_mode
         read(unit=lun,fmt=*) conf %user_modes% dcomp_mode
         read(unit=lun,fmt=*) conf % user_modes%acha_mode
         read(unit=lun,fmt=*) conf % user_modes%nlcomp_mode
         read(unit=lun,fmt=*) int_dummy
         conf % algo_group% run_sst = int_dummy == 1
         read(unit=lun,fmt=*) int_dummy
         conf % algo_group% run_cld = int_dummy == 1
         read(unit=lun,fmt=*) int_dummy
          conf % algo_group% run_aot= int_dummy == 1
         read(unit=lun,fmt=*) int_dummy
          conf % algo_group% run_erb= int_dummy == 1
         read(unit=lun,fmt=*)  int_dummy
         conf % algo_group% run_ash= int_dummy == 1
         read(unit=lun,fmt=*) conf % nwp_mode
         read(unit=lun,fmt=*)  int_dummy
         conf % do_smooth_nwp = int_dummy == 1
         read(unit=lun,fmt=*) conf % rtm_mode
         read(unit=lun,fmt=*)  int_dummy
         conf % do_process_cloudfree = int_dummy == 1
         read(unit=lun,fmt=*)  int_dummy
         conf % do_compress_out = int_dummy == 1
         read(unit=lun,fmt=*)  int_dummy
         conf % mask % read_aux = int_dummy == 1
         read(unit=lun,fmt=*)  int_dummy
         conf % mask % use_modis_clear = int_dummy == 1 
         read(unit=lun,fmt=*)  int_dummy
         conf % mask % use_prob_clear_res = int_dummy == 1 
         read(unit=lun,fmt=*)  int_dummy
         conf % mask % use_lrc = int_dummy == 1 
         read(unit=lun,fmt="(a)")   conf % mask %bayesian_mask_classifier
      end if
      close ( unit = lun )
      
   end subroutine read_user_configuration
   
   !
   !
   !
   subroutine print_it (self )
      class( conf_user_opt_type ) :: self
      integer :: i_file
      print*
      print*,'++++++++++++++++++++++++++++++++++++++++++ =++++++++++++++++++++++++++++++++++++++++'
      print*
      print*,'         configuration                 '
      print*
      print* , 'Ancil_path: ', trim(self % ancil_path)
      print* , 'Temp_path: ',trim(self % temp_path) 
      print*
      print* , 'User Algorithms Mode: '
      print* , ' MASK_MODE: ', self %user_modes% mask_mode
      print* , ' ACHA MODE: ', self %user_modes% acha_mode
      print* , ' DCOMP MODE: ', self %user_modes% dcomp_mode
      print*
      print* , 'Number of files: ', self % n_files
      do i_file =1 , self % n_files
         print*,'     File ',i_file,'    ',trim(self % file % infile (i_file)), ' is a ',trim(self % file % sensor_name (i_file)),' file'         
      end do
      print*
      print*,'++++++++++++++++++++++++++++++++++++++++++ =++++++++++++++++++++++++++++++++++++++++'
      print*
   
   
   end subroutine print_it

end module config_mod

