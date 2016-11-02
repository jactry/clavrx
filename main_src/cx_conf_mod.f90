! $Id$
module cx_conf_mod
   
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
      logical :: run_sasrab
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
      logical :: use_hres_sfctype
      logical :: use_hres_land_mask
      logical :: use_coast_mask
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
   
   type conf_output_type
      logical :: level2
      logical :: rtm
   end type conf_output_type
   
   type conf_file_type
      character (len=256) :: l1b_path
      character (len=256) :: out_path
      character (len=256) , allocatable, dimension (:)  :: infile
      character (len=256) , allocatable, dimension (:)  :: infile_fullpath
      integer , allocatable, dimension(:) :: wmo_id
      integer, allocatable, dimension (:)  :: ETsensor 
      integer, allocatable, dimension (:)  :: ETresolu 
   end  type conf_file_type
   
   type conf_algo_modes_type
      integer :: dcomp_mode
      integer :: mask_mode
      integer :: acha_mode
      integer :: nlcomp_mode
   end type conf_algo_modes_type
     
   type conf_main_type
      integer :: n_files
      integer :: expert_mode
      character (len=256) :: ancil_path
      character (len=256) :: temp_path
      type (conf_algo_modes_type ) :: user_modes
      type (conf_algo_modes_type ) :: updated_modes
      integer :: rtm_mode
      integer :: nwp_mode
      integer :: nav_mode
      logical :: do_smooth_nwp
      logical :: do_process_cloudfree
      logical :: do_compress_out
      integer :: num_scans_per_segment
      type (conf_process_type ) :: algo_group
      type (conf_diag_type ) :: diag
      type ( conf_sfc_type ) :: sfc
      type ( conf_mask_type) :: mask
      type ( conv_avhhr_type ) :: avhrr
      type ( conf_channel_type ) , dimension(48) :: chan
      type ( conf_output_type) :: out
      logical :: use_gzip     
      logical :: subset_pixel_hdf
      real, dimension(9) :: min_max
      
      type ( conf_file_type ):: file
      
      contains
      procedure, public :: set_config
      procedure, private :: read_user_configuration
      procedure, private :: read_input_files
      procedure, public :: update_modes
   end type conf_main_type
   
   enum , bind(C)
      enumerator :: ETsensor_invalid = 0
      enumerator :: ETsensor_avhrr  =1
      enumerator :: ETsensor_aqua 
      enumerator :: ETsensor_terra
      enumerator :: ETsensor_seviri 
      enumerator :: ETsensor_goes_first =5
      enumerator :: ETsensor_goes_mop =5
      enumerator :: ETsensor_goes_il =6
      enumerator :: ETsensor_goes_last =6
      enumerator :: ETsensor_viirs 
      enumerator :: ETsensor_mtsat
      enumerator :: ETsensor_aqua_mac 
      enumerator :: ETsensor_iff_viirs
      enumerator :: ETsensor_iff_aqua
      enumerator :: ETsensor_ahi
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
   
      use file_tools, only: get_lun &
         , file_nr_lines, file_test
      
      
      class ( conf_main_type ) :: conf
      character ( len =256 ) :: input_file = 'file_list'
      integer :: n_files
      integer :: i_file
      integer :: ios
      integer :: ppp
      
      if ( .not. file_test(trim(input_file))) return
     
      ppp = file_nr_lines (trim(input_file))
      
      n_files = ppp - 2
   
      conf % n_files = n_files
      allocate ( conf % file % infile (n_files) ) 
      allocate ( conf % file % infile_fullpath (n_files) )
      allocate ( conf % file % ETsensor ( n_files ) )
      allocate ( conf % file % wmo_id ( n_files ) )
      
      
       lun = get_lun()
      open(unit=lun,file = trim(input_file ),status="old",action="read",iostat=ios)
      read(unit=lun,fmt="(a)")  conf % file % l1b_path
      read(unit=lun,fmt="(a)")  conf % file % out_path
     
      do i_file =1 , n_files
         read(unit=lun,fmt="(a)")  conf % file % infile (i_file)
         
         conf % file % ETsensor (i_file) &
                 &  =  sensor_from_filename (trim(conf % file % infile (i_file)) , conf % file % l1b_path )
         conf % file % infile_fullpath (i_file) = trim(conf % file % l1b_path)//trim(conf % file % infile (i_file))
        
      end do
      
      close ( unit = lun )
      
      
   end subroutine read_input_files 
   
   !
   !  determine sensor 
   !
   function sensor_from_filename (file , path ) result (sensor)  
      use cx_ssec_areafile_mod
      character (len =* ) , intent(in) :: file
      character (len =* ) , intent(in) :: path
      integer  :: sensor
      integer :: resolu
      integer :: sat_id_num
           
      sensor = ETsensor_invalid
      resolu = ETresolu_invalid
      
      !-  test for modis
      if ( index( file , 'MYD') /= 0) sensor = ETsensor_aqua
      if ( index( file , 'MOD') /= 0) sensor = ETsensor_terra
      if ( index( file , 'MAC02')  /= 0) sensor = ETsensor_aqua_mac
      if ( index( file , 'HS_H08') /= 0) sensor = ETsensor_ahi
      
      if ( index( file , 'D02SSH') /= 0 ) resolu = ETresolu_5km
      if ( index( file , 'D021KM') /= 0 ) resolu = ETresolu_1km
      
      if ( sensor /= ETsensor_invalid) return
      
      !- test for VIIRS
      if ( index( file , 'GMTCO')  /= 0) sensor = ETsensor_viirs
      if ( sensor /= ETsensor_invalid) return
      
      ! - test for PEATE IFF VIIRS/MODIS data set
      if ( index( file , 'IFF_npp')  /= 0) sensor = ETsensor_iff_viirs
      if ( index( file , 'IFF_aqua') /= 0 ) sensor = ETsensor_iff_aqua
      
      ! - test if this is an area file
      if (   is_area_file ( trim(path)//file , sat_id_num ) ) then
         
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
      
      if ( sensor == ETsensor_invalid ) sensor = ETsensor_avhrr
     
   end function sensor_from_filename
   
   ! called from process_clavrx
   !
   !
   subroutine set_config ( conf)
      class ( conf_main_type ) :: conf
      
      call conf % read_user_configuration ()
      call conf % read_input_files
     
      
   end subroutine set_config 
   
   !
   !   copy user options to clavrx options
   !  and check if these are possible for this 
   !  sensor and file
   !
   subroutine update_modes ( this, i_file )
      class ( conf_main_type ) :: this
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
         print*, 'somehow wrong'
      
      end select 
      
   end subroutine update_modes
   
  !
  !
  !
   subroutine read_user_configuration ( conf )
      
      use cx_read_ascii_line_mod
      
      class ( conf_main_type ) :: conf
      character ( len =256 ) :: option_file = 'clavrx_options'
      integer :: lun
      integer :: ios0
      integer :: errstat
      integer :: int_dummy
      logical :: dum_l6(6)
      integer :: i
    
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
         read(unit = lun,fmt="(a)") conf % ancil_path
         read(unit = lun,fmt="(a)") conf % temp_path
         read(unit = Lun,fmt=*) conf % expert_mode
         
         if ( conf % expert_mode  == 0 )  then
            close(unit = Lun)
            return
         end if
                 
         read(unit=lun,fmt=*) conf %user_modes% mask_mode
         read(unit=lun,fmt=*) conf %user_modes% dcomp_mode
         read(unit=lun,fmt=*) conf % user_modes%acha_mode
         read(unit=lun,fmt=*) conf % user_modes%nlcomp_mode
         
         if ( conf % expert_mode  <= 1 )  then
            close(unit= Lun)
            return
         end if
         
         call read_ascii_line ( lun,  conf % out % level2)
         call read_ascii_line ( lun,  conf % out % rtm )
         call read_ascii_line ( lun,  conf % algo_group % run_cld )
         
         call read_ascii_line ( lun , conf % num_scans_per_segment)
         call read_ascii_line ( lun , conf % algo_group % run_sasrab )
         call read_ascii_line ( lun , conf % nwp_mode )
         call read_ascii_line ( lun , conf % rtm_mode )
         call read_ascii_line ( lun , conf % nav_mode )
         call read_ascii_line ( lun , conf % do_compress_out )
         
         call read_ascii_line ( lun , conf % mask % read_aux )
         call read_ascii_line ( lun , conf % mask % bayesian_mask_classifier )
         
         if ( conf % expert_mode  <= 2 )  then
            close(unit= Lun)
            return
         end if
         
         
         call read_ascii_line ( lun , conf % sfc % use_seebor )
         call read_ascii_line ( lun , conf % sfc % use_hres_sfctype )
         
         call read_ascii_line ( lun , conf % sfc % use_hres_land_mask)
         call read_ascii_line ( lun , conf % sfc % use_coast_mask)
         
         call read_ascii_line ( lun , conf % sfc % use_hres_elev)
         
         call read_ascii_line ( lun , conf % sfc % use_volc_mask)
         call read_ascii_line ( lun , conf % sfc % snow_mask_mode)
         call read_ascii_line ( lun , conf % sfc % use_dark_compo)
         
          if ( conf % expert_mode  <= 3 )  then
            close(unit= Lun)
            return
         end if
         
         call read_ascii_line ( lun , conf % avhrr % active_rfl_cal_1b )
         call read_ascii_line ( lun , conf % avhrr % use_therm_cal_1b )
         
         if ( conf % expert_mode  <= 4 )  then
            close(unit= Lun)
            return
         end if
         
         call read_ascii_line ( lun , conf % mask % use_lrc )
         call read_ascii_line ( lun , conf % do_smooth_nwp )
         call read_ascii_line ( lun , conf % do_process_cloudfree )
         call read_ascii_line ( lun, conf % min_max )
         if ( conf % expert_mode  <= 5 )  then
            close(unit= Lun)
            return
         end if
         
         do i = 0,7      
            call read_ascii_line ( lun , dum_l6 )
            
            conf % chan ( (i * 6) + 1 : (i+1)* 6 ) % is_on = dum_l6
         end do
         
      end if
      close ( unit = lun )
      
      
    
      
   end subroutine read_user_configuration

end module cx_conf_mod


