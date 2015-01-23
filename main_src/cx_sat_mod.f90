! $Header: https://svn.ssec.wisc.edu/repos/aw_clavrx/trunk/sat_mod.f90 18 2014-10-23 20:40:18Z awalther $
!
module cx_sat_mod
  ! use sensor_mod
   use date_tools_mod
   
   
   implicit none
   
   real , parameter :: PI = 3.14159265359
   real , parameter :: DTOR = PI / 180.
   
   integer , dimension(16) , private , parameter :: & 
      modis_chn_list = [ 8 , 9 , 10 , 12 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
   
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
      character (len = 255) ::  Ancil_data_dir
      integer :: ETsensor
      integer :: num_pix
      integer :: num_elem
      integer :: num_seg
      integer :: num_elem_per_seg_general
      integer :: num_elem_this_seg
      integer :: idx_seg
      logical , dimension(42) :: chan_on
      type ( sat_data_type) , dimension(42) :: chn
      type ( geo_data_type) :: geo
      type ( date_type ) :: start_time
      type ( date_type ) :: end_time
   !   type ( sensor_type) :: sensor_const
      
   contains
      procedure :: init   
      procedure :: set_dimensions 
      procedure :: read_l1b
      procedure :: deallocate_all
      procedure :: determine_channels
      
   end type sat_main_type
   
   
    

contains

   subroutine init ( this , conf , i_file) 
      use config_mod
      class ( sat_main_type) :: this
      type ( conf_user_opt_type ) , intent ( in ) :: conf 
      integer , intent (in) :: i_file
     
      this % file = conf % file % infile (i_file)
      this % ETsensor = conf % file % ETsensor (i_file)
      this % l1b_path = conf % file % l1b_path
      this % ancil_data_dir =  conf % ancil_path
      
    !  associate ( sensor => this % sensor_const)
      !   call sensor % populate(this%file, this%l1b_path)
    ! end associate
      
    !  call this % set_dimensions
      
    !  print*,'   ==   Sensor is ===>   ',this%sensor_const%sensor_name
      
   end subroutine init 
   !
   !
   !
   subroutine set_dimensions ( this )
      use config_mod
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
          
               this % num_pix = 3200
               this % num_elem = 768
               
               this % num_seg = 4
               this % num_elem_per_seg_general = 200  ! - later via option file
              
             !  call time_0 % print_data  
             !  call time_1 % print_data 
                
         end select
      end associate
   
   end subroutine set_dimensions
   
   ! ===============================================
   !
   ! ===============================================
   
   subroutine read_l1b ( this , idx_seg,areastr,navstr,Nrec_Avhrr_Header,time_since_launch)
      use pixel_common
      use sensor_module
     
           
      class ( sat_main_type) :: this
      integer , intent(in) :: idx_seg
      TYPE (AREA_STRUCT) :: AREAstr
      TYPE (GVAR_NAV)    :: NAVstr
      integer :: ierror_level1b
      integer(kind=int4):: Nrec_Avhrr_Header
      real(kind=real4):: Time_Since_Launch
      
      integer :: ii
      
      print*,'====> ',idx_seg
      ! - this is old style and has to be rewritten
      call READ_LEVEL1B_DATA(Image%Level1b_Full_Name,idx_seg, &
                                Time_Since_Launch,AREAstr,NAVstr,Nrec_Avhrr_Header,Ierror_Level1b)
                                
         
         
         if (Ierror_Level1b /= 0) then
            print *, EXE_PROMPT, "ERROR:  Error reading level1b, skipping this file"
            
         end if
         
      this % file = 'the file'
      this % num_pix = Image%Number_Of_Lines
      this % num_elem_this_seg = Image%Number_Of_Lines_Per_Segment
      
      do ii = 1, 42
         
         if ( allocated ( ch(ii) % rad_toa)) this % chn(ii) % rad = ch(ii) % rad_toa
         if ( allocated ( ch(ii) % ref_toa)) this % chn(ii) % ref = ch(ii) % ref_toa
         if ( allocated ( ch(ii) % bt_toa)) this % chn(ii) % bt = ch(ii) % bt_toa
         
      end do   
      if ( .not. allocated  (nav%lat_1b) ) then
			print*,'nav lon is to existing'
			stop
		end if
		
		
      if ( allocated (nav%lat_1b)) this %  geo % lat = nav % lat_1b
      if ( allocated (nav%lon_1b)) this %  geo % lon = nav % lon_1b
      if ( allocated (geo % solzen)) this %  geo % sol_zen = geo % solzen
      if ( allocated (geo % satzen)) this %  geo % sat_zen = geo % satzen
      if ( allocated (geo % relaz)) this %  geo % rel_az = geo % relaz
      if ( allocated (geo % scatangle)) this %  geo % scat_zen = geo % scatangle
      if ( allocated (geo % glintzen)) this %  geo % glint_zen = geo % glintzen
      !if ( allocated (nav%lat_1b)) this %  geo % lat = nav % airmass
   
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
     
      
      do i_mband = 1 , 42 
        
         if ( allocated ( this % chn (i_mband)  % ref ) ) deallocate ( this % chn (i_mband)  % ref )
         if ( allocated ( this % chn (i_mband)  % rad ) ) deallocate ( this % chn (i_mband)  % rad )
         if ( allocated ( this % chn (i_mband)  % bt ) ) deallocate ( this % chn (i_mband)  % bt )
         
      end do
   
   end subroutine deallocate_all
   
   !
   !
   !
   subroutine determine_channels  ( this, &
             dcomp_mode, acha_mode , mask_mode )
   
   !   use sensor_mod   
      
      class ( sat_main_type) :: this
      
      
      integer :: dcomp_mode
      integer :: acha_mode
      integer :: mask_mode
      

      this % chan_on(:) = .false.
      
 !     select case ( this % ETsensor)
   
 !     case ( ETsensor_avhrr )
 !        this % chan_on (1) = .true.
 !        this % chan_on (2) = .true.
  !       this % chan_on (20) = .true.
 !        this % chan_on (31) = .true.
 !        this % chan_on (32) = .true.
   
 !     case ( ETsensor_viirs)
 !        this % chan_on (1) = .true.
 !        this % chan_on (2) = .true.
  !       this % chan_on (6) = .true.
  !       this % chan_on (20) = .true.
  !       this % chan_on (26) = .true.
  !       this % chan_on (29) = .true.
  !       this % chan_on (31) = .true.
  !       this % chan_on (32) = .true.
  !    end select
   end subroutine determine_channels
   


end module cx_sat_mod
