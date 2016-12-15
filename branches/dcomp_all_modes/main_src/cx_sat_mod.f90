! $Header$
!
module cx_sat_mod
   use sensor_mod
   use date_tools_mod, only: &
      date_type
   
   
   implicit none
   
   real , parameter :: PI = 3.14159265359
   real , parameter :: DTOR = PI / 180.
   
   integer , dimension(16) , private , parameter :: & 
      modis_chn_list = [ 8 , 9 , 10 , 12 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
    
    integer , dimension(16) , private , parameter :: &
    modis_chn_list_ahi = [ 3 , 4 , 1 , 2 , 6 , 7 , 20 , 43 ,  27 , 28 , 29 , 30 , 44 , 31 , 32 , 33 ]
   type sat_data_type
      real, allocatable, dimension(:,:) :: rad
      real, allocatable, dimension(:,:) :: ref
      real, allocatable, dimension(:,:) :: bt
   end type sat_data_type
   
   type geo_data_type
      real, allocatable, dimension(:,:) :: lat
      real, allocatable, dimension(:,:) :: lon
      real, allocatable, dimension(:,:) :: sol_zen
      real, allocatable, dimension(:,:) :: sat_zen
      real, allocatable, dimension(:,:) :: rel_az
      real, allocatable, dimension(:,:) :: scat_zen
      real, allocatable, dimension(:,:) :: glint_zen
      real, allocatable, dimension(:,:) :: airmass_sol
   end  type geo_data_type
   
   type sat_main_type
      character (len = 255) ::  file
      character (len = 255) ::  l1b_path
      character ( len =1000) :: path_file
      character (len = 255) ::  Ancil_data_dir
      integer :: ETsensor
      integer :: num_pix
      integer :: num_elem
      integer :: num_seg
      integer :: num_elem_per_seg_general
      integer :: num_elem_this_seg
      integer :: idx_seg
      logical , dimension(45) :: chan_on
      type ( sat_data_type) , dimension(45) :: chn
      type ( geo_data_type) :: geo
      type ( date_type ) :: start_time
      type ( date_type ) :: end_time
      type ( sensor_type) :: sensor_const
      
   contains
      procedure :: init   
      procedure :: set_dimensions 
      procedure :: read_l1b
      procedure :: deallocate_all
      
      
   end type sat_main_type
   
   
    

contains

   function init ( this , conf , i_file ) result (ok) 
      use cx_conf_mod, only: &
         conf_main_type
      use file_tools
      class ( sat_main_type) :: this
      type ( conf_main_type ) , intent ( in ) :: conf 
      integer , intent (in) :: i_file
      logical   :: ok
      integer :: i
     
      ok = .false.
      
      this % file = conf % file % infile (i_file)
      
      this % ETsensor = conf % file % ETsensor (i_file)
      this % l1b_path = conf % file % l1b_path
      this % ancil_data_dir =  conf % ancil_path
      this % path_file = trim(this % l1b_path)//trim(this % file) 
     
      if ( file_test (trim(this % path_file)) ) then
            call this % sensor_const % populate(this%file, this%l1b_path)
            call this % set_dimensions
            ok = .true.
      end if   
      
      !AW-TODO has to come from config
      do i=1, 45
         this % chan_on(i) = conf % chan (i) %is_on
      end do
   end function init 
   !
   !
   !
   subroutine set_dimensions ( this )
      class ( sat_main_type) :: this
      integer :: year, month , day , start_hour, start_minute
      integer :: start_sec, end_hour , end_minute, end_sec, orbit
      character ( len =100)::infile
    
      
      associate ( time_0 => this %  start_time , &
                       time_1 => this % end_time     )
     
         select case ( this % ETsensor )
            case (ETsensor_aqua)
      
            case (ETsensor_viirs)
               ! --- Read data from the file name
               infile = this % file
         
               read(Infile(12:15), fmt="(I4)") year
               read(Infile(16:17), fmt="(I2)") month
               read(Infile(18:19), fmt="(I2)") day
               read(Infile(22:23), fmt="(I2)") start_hour
               read(Infile(24:25), fmt="(I2)") start_minute
               read(Infile(26:27), fmt="(I2)") start_sec
               read(Infile(31:32), fmt="(I2)") end_hour
               read(Infile(33:34), fmt="(I2)") end_minute
               read(Infile(35:36), fmt="(I2)") end_sec
               read(Infile(40:44), fmt="(I5)") orbit
         
            
          
               this % num_pix = 3200
               this % num_elem = 768
               
               this % num_seg = 4
               this % num_elem_per_seg_general = 200  ! - later via option file
              
            case ( ETsensor_ahi )
               infile = this % file
               print*,infile
               
               read(Infile(8:11), fmt ="(I4)") year
               read(Infile(12:13), fmt="(I2)") month
               read(Infile(14:15), fmt="(I2)") day
               read(Infile(17:28), fmt="(I2)") start_hour
               read(Infile(19:20), fmt="(I2)") start_minute
               
               
               start_sec = 0
               end_hour = start_hour
               end_minute = start_minute+10
               end_sec = 0
               
               this % num_pix = 5500
               this % num_elem = 5500
               
               this % num_seg = 28
               this % num_elem_per_seg_general = 200  ! - later via option file
               
                
         end select
         print*,'gggg==> ', this % ETsensor
         print*,year,month,day,start_hour,start_minute,start_sec
         
         call time_0 % set_date (  &
                  &  year = year &
                  &  , month = month &
                  & , day = day  &
                  & , hour = start_hour &
                  & , minute = start_minute &
                  & , second = start_sec )
                  
          call time_1 % set_date (  &
                  &  year = year &
                  &  , month = month &
                  & , day = day  &
                  & , hour = end_hour &
                  & , minute = end_minute &
                  & , second = end_sec )    
         
         
      end associate
   
   end subroutine set_dimensions
   
   ! ===============================================
   !
   ! ===============================================
   
   subroutine read_l1b ( this , idx_seg , error_status)
   
      use viirs_read_mod , only : &
             viirs_data_config &
           , viirs_data_out &
           , get_viirs_data
      use cx_read_ahi_mod, only : &
         ahi_config_type &
         , ahi_data_out_type &
         , get_ahi_data
           
  !    use planck_mod , only : &
    !     planck_rad2bt     
         
      use viewing_geometry_module, only : &
         glint_angle  , scattering_angle 
           
      class ( sat_main_type) :: this
      integer , intent(in) :: idx_seg
      logical , intent(out) :: error_status 
      
      type ( ahi_config_type ) :: ahi_c
      type ( ahi_data_out_type ) :: ahi_data
      type ( viirs_data_config )  :: v_conf
      type ( viirs_data_out )  :: out
      integer :: i_mband
      
      error_status = .false.
      
     
      if ( this % ETsensor /= ETsensor_viirs .and. this % ETsensor /= ETsensor_ahi ) then
           print*,'non-viirs and non-ahi sensors have to be impolemented step-by-step '
           error_status = .true.
         return
      end if
      
      ! - configure viirs interface
      this % num_elem_this_seg = this % num_elem_per_seg_general
      if ( idx_seg == this % num_seg) then
         this % num_elem_this_seg = this % num_elem - (idx_seg - 1 ) * this % num_elem_per_seg_general    
      end if
      
      allocate ( this % geo  % lat        (this % num_pix, this % num_elem_this_seg))
      allocate ( this % geo  % lon        (this % num_pix, this % num_elem_this_seg))
      allocate ( this % geo  % sol_zen    (this % num_pix, this % num_elem_this_seg))
      allocate ( this % geo  % sat_zen    (this % num_pix, this % num_elem_this_seg))
      allocate ( this % geo  % scat_zen   (this % num_pix, this % num_elem_this_seg))
      allocate ( this % geo  % rel_az     (this % num_pix, this % num_elem_this_seg))
      allocate ( this % geo  % glint_zen  (this % num_pix, this % num_elem_this_seg))
      allocate ( this % geo  % airmass_sol(this % num_pix, this % num_elem_this_seg))
      
      
     
      select case ( this % ETsensor )
      
      case (ETsensor_viirs)
      
         v_conf % chan_on_rfl_mband = this % chan_on (modis_chn_list)        
         v_conf % file_gmtco_base =  trim( this % file)
         v_conf % dir_1b = trim(this %l1b_path )      
         v_conf % offset (1) = 1
         v_conf % offset (2)  = (idx_seg -1 ) * this % num_elem_per_seg_general + 1
         v_conf % count (1) = 3200
         v_conf % count (2)  = this % num_elem_this_seg
         v_conf % ancil_data_dir = this % ancil_data_dir
     
         call get_viirs_data ( v_conf, out )
         
         this % geo % lat =      out % geo % lat
         this % geo % lon =      out % geo % lon
         this % geo % sol_zen =  out % geo % solzen
         this % geo % sat_zen =  out % geo % satzen
         this % geo % rel_az =   out % geo % relaz
         this % geo % glint_zen =  glint_angle ( this % geo % sol_zen , this % geo % sat_zen , this % geo % rel_az )
         this % geo % airmass_sol = 1/cos (this % geo % sat_zen * DTOR) + 1/ cos (this % geo % sol_zen * DTOR) 
         this % geo % scat_zen = scattering_angle ( this % geo %sol_zen, this % geo %sat_zen, this % geo %rel_az ) 
         
         do i_mband = 1, 16
           
            if ( .not. this % chan_on (modis_chn_list(i_mband))  ) cycle
                
            if ( i_mband < 12 ) then
               allocate (     this % chn (modis_chn_list(i_mband))  % ref (3200,this % num_elem_this_seg))
               this % chn (modis_chn_list(i_mband))  % ref =  out % mband (i_mband) % ref   
            else
               
               allocate ( this % chn (modis_chn_list(i_mband))  % rad (3200,this % num_elem_this_seg))
               allocate ( this % chn ( modis_chn_list(i_mband)) % bt (3200, this % num_elem_this_seg))
               
               this % chn (modis_chn_list(i_mband))  % rad =  out % mband (i_mband) % rad 
        !    call planck_rad2bt (this % chn ( modis_chn_list(i_mband)) % rad, 'VIIRS' &
          !     , modis_chn_list(i_mband) , this % chn (modis_chn_list(i_mband) )% bt  )  
            end if
         
         end do
      
      case ( ETsensor_ahi)
         ahi_c % file_base = trim(this % file)
         ahi_c % data_path = trim ( this %l1b_path )
         ahi_c % h5_offset = [0,(idx_seg -1 ) * this % num_elem_per_seg_general ]
         ahi_c % h5_count = [5500,this % num_elem_this_seg]    
         ahi_c % chan_on (:) = .true.
         
         call get_ahi_data ( ahi_c , ahi_data )
         this % geo % lat =      ahi_data % geo % lat
         this % geo % lon =      ahi_data % geo % lon

         do i_mband = 1,16 

            if ( .not. this % chan_on (modis_chn_list_ahi(i_mband))  ) cycle
            if ( i_mband < 7 ) then
               allocate ( this % chn (modis_chn_list_ahi(i_mband))  % ref (this % num_pix,this % num_elem_this_seg))
               this % chn (modis_chn_list_ahi(i_mband))  % ref = ahi_data % chn ( i_mband) % ref
            else
               
               allocate ( this % chn (modis_chn_list_ahi(i_mband))  % rad (this % num_pix,this % num_elem_this_seg))
               this % chn (modis_chn_list_ahi(i_mband))  % rad = ahi_data % chn ( i_mband) % rad
            end if   
         end do
         
         call ahi_data % deallocate_all
         
         
      end select
      
   
     
    
    
     
   
   
   end subroutine
   
   ! ======================================
   !
   ! ======================================   
   subroutine deallocate_all ( this )
      class ( sat_main_type) :: this
      integer :: i_mband
      
      if ( allocated ( this % geo  % lon ) ) deallocate ( this % geo  % lat )
      if ( allocated ( this % geo  % lon ) ) deallocate ( this % geo  % lon )
      if ( allocated ( this % geo  % sol_zen ) ) deallocate ( this % geo  % sol_zen )
      if ( allocated ( this % geo  % sat_zen ) ) deallocate ( this % geo  % sat_zen )
      if ( allocated ( this % geo  % scat_zen ) ) deallocate ( this % geo  % scat_zen )
      if ( allocated ( this % geo  % rel_az ) ) deallocate ( this % geo  % rel_az )
      if ( allocated ( this % geo  % glint_zen ) ) deallocate ( this % geo  % glint_zen )
      if ( allocated ( this % geo  % airmass_sol ) ) deallocate ( this % geo  % airmass_sol )
     
      
      do i_mband = 1 , 44
       
         if ( allocated ( this % chn (i_mband)  % ref ) ) deallocate ( this % chn (i_mband)  % ref )
         if ( allocated ( this % chn (i_mband)  % rad ) ) deallocate ( this % chn (i_mband)  % rad )
         if ( allocated ( this % chn (i_mband)  % bt ) ) deallocate ( this % chn (i_mband)  % bt )
         
      end do
   
   end subroutine deallocate_all
   

   


end module cx_sat_mod
