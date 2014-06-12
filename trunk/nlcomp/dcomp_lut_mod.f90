! $Id: dcomp_lut_mod.f90 422 2014-06-12 16:13:30Z awalther $
!
! HISTORY: 06/05/2014: change of filename 
!
module dcomp_lut_mod
   
   use dcomp_math_tools_mod, only:&
      dcomp_interpolation_weight &
      , interpolate_2d
      
	use file_tools, only: &
      file_test
   
   
    use ica_f90_hdf_sds , only : &
         & hdf_sds &
         , hdf_data &
         , hdf_get_file_sds &
         , MAXNCNAM
   
   implicit none
 
   private
   
   integer, parameter :: NUM_PHASE = 2
   integer, parameter :: NUM_CHN = 43
   
   logical :: is_initialized = .false.
   character(10) :: sensor_initialized
   
   
   type lut_data_type
      logical :: is_set
      logical :: has_sol
      logical :: has_ems
      character (len = 300 ) :: file
      character (len = 300 ) :: file_ems
      real, allocatable :: cld_sph_alb ( : , : )
      real, allocatable :: cld_trn ( : , : , : )
      real, allocatable :: cld_alb ( : , : , : )
      real, allocatable :: cld_refl( : , : , : , : , : )
      real, allocatable :: cld_ems ( : , : , : )
      real, allocatable :: cld_trn_ems ( : , : , : )
      
      contains
      procedure :: read_hdf => lut_data__read_hdf
      procedure :: alloc => lut_data__alloc
      procedure :: dealloc => lut_data__dealloc
      
   end type lut_data_type
   
   type lut_chn_type
      logical :: is_set
     
      type ( lut_data_type ) , dimension ( NUM_PHASE ) :: phase
   end type lut_chn_type
   
   type lut_dim_type
      logical :: is_set
      real , allocatable :: sat_zen ( : )
      real , allocatable :: sol_zen ( : )
      real , allocatable :: rel_azi ( : )
      real , allocatable :: cod ( : )
      real , allocatable :: cps ( : )  
      integer :: n_sat_zen
      integer :: n_sol_zen
      integer :: n_rel_azi
      integer :: n_cod
      integer :: n_cps
      
      contains
      procedure :: alloc => lut_dim__alloc
      procedure :: dealloc => lut_dim__dealloc
       
   end type lut_dim_type
   
   type :: lut_type
      logical :: is_set
      type ( lut_chn_type ) , dimension ( NUM_CHN ) :: channel 
      type ( lut_dim_type ) ::  dims
      character ( len = 300) :: lut_path
      character ( len = 20 ) :: sensor = 'not_set'
      integer :: pos_sat, pos_sol, pos_azi
      
      contains
      procedure :: initialize => lut__initialize
      procedure :: set_angles => lut__set_angles
      procedure :: getProperty => lut__getProperty
      procedure :: get_data => lut__get_data
      procedure :: thick_cloud_rfl => lut_data__thick_cloud_rfl
      procedure :: init_dims => lut__init_dims
      procedure , private :: set_filename => lut__set_filename
      procedure :: clear_lut_memory => lut__clear_lut_memory
      
   end type lut_type
   
   ! - only one channel can be stored 
   type ( lut_type ) , public ,save , target :: lut_obj
   
   
   type , public :: lut_output
   
      real :: refl
      real :: trn_sol
      real :: trn_sat
      real :: albsph
      real :: alb
      real :: ems
      real :: trn_ems
      
      real :: dRefl_dcps      , dRefl_dcod
      real :: dtrans_sol_dcps , dtrans_sol_dcod
      real :: dTrans_sat_dcod , dTrans_sat_dcps
      real :: dsph_alb_dcod   , dSph_alb_dcps
      real :: dalb_dcod   , dalb_dcps
      real :: dEms_dcps      , dEms_dcod
      real :: dtrnEms_dcps      , dtrnEms_dcod
   
   
   end type lut_output
   

   real , save :: wgt_sat
   real , save :: wgt_sol
   real , save :: wgt_azi
   
   real, save :: sat_m
   real, save :: sol_m
   real, save :: azi_m
      
contains
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__set_filename ( self)
      class ( lut_type ) :: self
      character ( len = 3 ) , dimension(30) :: chan_string ='no'
      logical , dimension ( 43 ) :: has_ems_table = .false.
	   logical , dimension ( 43 ) :: has_sol_table = .false.
      
      integer :: i_chn , i_phase , i 
      character ( len = 3 ) , dimension(2)   :: phase_string = [ 'wat',  'ice' ]
      integer :: n_channels = 43
      character ( len =200) :: sensor_identifier
      
      ! mapping sensor channel emis yes/no
      
      sensor_identifier = trim(self % lut_path) & 
                           & //  trim ( self % sensor )
     
	   sensor_block: select case ( trim(self % sensor))
	   case ('Meteosat-8','Meteosat-9','Meteosat-10') sensor_block
		   has_sol_table(1:2) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(2) = '2'
			chan_string(6) = '3'
			chan_string(20) = '4'
      
      case ('NOAA-05','NOAA-06','NOAA-07','NOAA-08','NOAA-09', 'NOAA-10','NOAA-11','NOAA-12', 'NOAA-14')  sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(20) = '3b'
        
		case ('NOAA-15','NOAA-16', 'NOAA-17','NOAA-18','NOAA-19','METOP-A','METOP-B')  sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(6) = '3a'
			chan_string(20) = '3b'		
      
      case ('GOES-08','GOES-09','GOES-10', 'GOES-11' , 'GOES-12' , 'GOES-13',  'GOES-14', 'GOES-15','COMS-1'  )   sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(20) = '2'
      
      case ('MODIS-AQUA', 'MODIS-TERRA')    sensor_block
         has_sol_table(1:2) = .true.
			has_sol_table(5:7) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(2) = '2'
			chan_string(5) = '5'
			chan_string(6) = '6'
			chan_string(7) = '7'
			chan_string(20) = '20'
         sensor_identifier = trim(self % lut_path) //'MODIS' 
                          
         
      case('GOES-16')  sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '2'
			chan_string(6) = '5'
			chan_string(20) = '7'
         
      case('ABI') sensor_block
         has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '2'
			chan_string(6) = '5'
			chan_string(20) = '7'
            
      case ('AATSR')   sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(6) = '6'
			chan_string(20) = '20'
             
      case ('VIIRS')   sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(5) = .true.
			has_sol_table(6) = .true.
			has_sol_table(7) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '5'
			chan_string(5) = '8'
			chan_string(6) = '10'
			chan_string(7) = '11'
			chan_string(20) ='12'
         
      case ('MTSAT-1R')   sensor_block
         has_sol_table(1) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(20) = '5'  
         
      case ('MTSAT-2')   sensor_block
         has_sol_table(1) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(20) = '5'  	
                
      case default
          print*,'add sensor in dcomp_data_pool_mod.f90 routine populate...', trim(self%sensor)
          stop
	  
	   end select sensor_block
     
      loop_channel : do i_chn = 1 , n_channels
         if ( .not. has_sol_table ( i_chn ) )  cycle
         loop_phase: do i_phase = 1 , 2
            self%channel(i_chn)%phase(i_phase)%file = trim(sensor_identifier) & 
                       
                       & // '_ch'//trim ( chan_string ( i_chn ) ) &
                       & //'_ref_lut_'//phase_string(i_phase)//'_cld.hdf'
                   
                  
            if ( has_ems_table(i_chn) ) then  
               self%channel(i_chn)%phase(i_phase)%file_ems = trim(sensor_identifier) & 
                     
                       & // '_ch'//trim ( chan_string ( i_chn ) ) &
                       & //'_ems_lut_'//phase_string(i_phase)//'_cld.hdf'
            end if         
         end do loop_phase
      end do loop_channel
      
      do i =1,2 
         self % channel(:) % phase(i) % has_ems = has_ems_table
	  end do
   
   end subroutine lut__set_filename
   
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__initialize ( self, sensor , ancil_path )
      class ( lut_type ) :: self
      character ( len = * ) , intent(in) :: sensor
      character ( len = * ) , intent(in), optional :: ancil_path
      character ( len =300) :: file
      character ( len =20) :: host
      
      ! - check if sensor is already initialized
      if ( self % sensor == sensor ) then
         return
      end if
      
      call getenv("HOST",host)
      print*,'initialized ', sensor
      
      ! - some lut paths
      self % lut_path = '/DATA/Ancil_Data/clavrx_ancil_data/luts/cld/'
		if ( host(1:4) == 'saga' ) self % lut_path = '/data/Ancil_Data/clavrx_ancil_data/luts/cld/' 
      if ( present(ancil_path)) self % lut_path = trim(ancil_path)
      self % sensor = trim(sensor)
      
      ! - set filenames
      call self % set_filename
      ! - clear memory for new sensor   
      call self % clear_lut_memory
      call self % init_dims( )
      
   end subroutine lut__initialize
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__set_angles ( self, sat , sol , azi )
      class ( lut_type ) :: self
      real , intent(in) , optional :: sat
      real , intent(in) , optional :: sol   
      real , intent(in) , optional :: azi
      
      if ( present(sat) ) sat_m = sat
      if ( present(sol) ) sol_m = sol
      if ( present(azi) ) azi_m = azi
      
      ! - compute pos and weights
      call dcomp_interpolation_weight(self%dims%n_sat_zen , sat , self%dims%sat_zen &                  
                   &, near_index = self % pos_sat  )
      call dcomp_interpolation_weight(self%dims%n_sol_zen , sol , self%dims%sol_zen &                  
                   &, near_index = self % pos_sol  )
      call dcomp_interpolation_weight(self%dims%n_rel_azi , azi , self%dims%rel_azi &                  
                   &, near_index = self % pos_azi  )  
                  
   
   end subroutine lut__set_angles
   
   ! --------------------------------------------------------------------------
   !
   !   clears memory of LUT_data object
   !   
   ! --------------------------------------------------------------------------   
   subroutine lut__clear_lut_memory ( self )
   class ( lut_type ) , target :: self
   integer :: idx_phase , idx_chn
   type (lut_data_type), pointer :: data_loc => null()
   
   
   do idx_phase =1, NUM_PHASE
      do idx_chn =1, NUM_CHN
         data_loc => self % channel ( idx_chn ) % phase ( idx_phase)
         call data_loc % dealloc()
         
      end do
   end do   
   
   
   
   end subroutine lut__clear_lut_memory
   
   !  --------------------------------------------------------------------------
   !   PURPOSE : return data from LUT
   !   input: channel, phase, cod (log10 ), cps(log10)
   !   output : transmission, reflectance, cloud albedo, spherical abedo
   !          several derivates
   !   derivates we need:
   !
   !       drefl_dcps
   !       drefl_dcod
   !       
   !
   !  -------------------------------------------------------------------------
   subroutine lut__get_data ( self, idx_chn , idx_phase , cod_log10, cps_log10 &
                              & , out )
                              
       class ( lut_type ) , target :: self
      integer , intent(in) :: idx_chn
      integer , intent(in) :: idx_phase
      real, intent(in) :: cod_log10
      real, intent(in) :: cps_log10
      
      type ( lut_output ) :: out
      
      integer :: pos_cod,pos_cps
      real :: wgt_cod, wgt_cps
      real , save :: cod_log10_saved
      
      
      real :: rfl_cld_2x2 (2,2)
      real :: trn_sol_cld_2x2 (2,2)
      real :: trn_sat_cld_2x2 (2,2)
      real :: albsph_cld_2x2 (2,2)
      real :: ems_cld_2x2 ( 2,2)
      real :: trn_ems_cld_2x2 ( 2,2)
      real :: alb_cld_2x2 ( 2,2)
      real :: ref_diff , cod_diff
      
      real :: dRefl_dcps      , dRefl_dcod
      real :: dtrans_sol_dcps , dtrans_sol_dcod
      real :: dTrans_sat_dcod , dTrans_sat_dcps
      real :: dsph_alb_dcod   , dSph_alb_dcps
      real :: dEms_dcps      , dEms_dcod
      
      
      type (lut_data_type), pointer :: data_loc => null()
      
      
      out % ems = -999.
       
       data_loc => self % channel ( idx_chn ) % phase ( idx_phase)
      
      ! test if this channel is available if not read it from hdf file     
      if ( .not. data_loc % is_set ) then 
         call data_loc % read_hdf
         data_loc % is_set  = .true.
      end if
      
      call dcomp_interpolation_weight(self%dims%n_cod, cod_log10,self%dims%cod &
         & , weight_out = wgt_cod, index_out= pos_cod)
      
      call dcomp_interpolation_weight(self%dims%n_cps, cps_log10,self%dims%cps &
         & , weight_out = wgt_cps, index_out= pos_cps)
         
         
      rfl_cld_2x2       = data_loc%cld_refl( pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sol,self%pos_sat,self%pos_azi)
      trn_sol_cld_2x2   = data_loc%cld_trn(pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sol)
      trn_sat_cld_2x2   = data_loc%cld_trn(pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sat)
      albsph_cld_2x2    = data_loc%cld_sph_alb(pos_cps:pos_cps+1,pos_cod:pos_cod+1)
      alb_cld_2x2       = data_loc%cld_alb(pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sol)
      
      ! - parameter for kernel computation  
      ref_diff = 0.2
	   cod_diff = 0.1
      
	   call interpolate_2d ( rfl_cld_2x2 , wgt_cps , wgt_cod , ref_diff , cod_diff , out % refl    &
         & , out % dRefl_dcps      , out % dRefl_dcod ) 
	   call interpolate_2d ( trn_sol_cld_2x2 , wgt_cps , wgt_cod , ref_diff , cod_diff , out % trn_sol  &
         & , out % dtrans_sol_dcps , out % dtrans_sol_dcod)
      call interpolate_2d ( trn_sat_cld_2x2 , wgt_cps , wgt_cod , ref_diff , cod_diff , out % trn_sat &
         & , out % dTrans_sat_dcod , out % dTrans_sat_dcps )
	   call interpolate_2d ( albsph_cld_2x2  , wgt_cps , wgt_cod , ref_diff , cod_diff , out % albsph &
         & , out % dsph_alb_dcod   , out % dSph_alb_dcps) 
      call interpolate_2d ( alb_cld_2x2  , wgt_cps , wgt_cod , ref_diff , cod_diff , out % alb &
         & , out % dalb_dcod   , out % dalb_dcps) 
      
      
      if ( data_loc % has_ems ) then
         
         ems_cld_2x2 = data_loc%cld_ems (pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sat)        
         call interpolate_2d ( ems_cld_2x2, wgt_cps , wgt_cod , ref_diff , cod_diff &
            , out % ems, out % dEms_dcps , out % dEms_dcod )
         
         trn_ems_cld_2x2 = data_loc%cld_trn_ems (pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sat)        
         call interpolate_2d ( trn_ems_cld_2x2, wgt_cps , wgt_cod , ref_diff , cod_diff &
            , out % trn_ems, out % dtrnEms_dcps, out % dtrnEms_dcod ) 
         
      end if  
      
      cod_log10_saved = cod_log10
   end subroutine lut__get_data
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__thick_cloud_rfl ( self , idx_chn , idx_phase, rfl, ems)
      class ( lut_type ) , target :: self
      integer , intent(in) :: idx_chn
      integer , intent(in) :: idx_phase
      
      real, intent(out) :: rfl(9)
      real, intent(out) :: ems(9)
      
      type (lut_data_type), pointer :: data_loc => null()
      
      data_loc => self % channel ( idx_chn ) % phase ( idx_phase)
      
      rfl = data_loc%cld_refl( :,29,self%pos_sol,self%pos_sat,self%pos_azi)
      if ( data_loc % has_ems ) then
         ems = data_loc%cld_ems (:,29,self%pos_sat)
      end if
   
   
   end subroutine lut_data__thick_cloud_rfl
   
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__read_hdf ( self )
      
      class ( lut_data_type ) :: self
      
      integer :: nsds
      character(len =MAXNCNAM), allocatable :: sds_name(:)
      type(hdf_sds) , allocatable, target :: sds(:)
      character(len =MAXNCNAM), allocatable :: sds_name_ems(:)
      type(hdf_sds) , pointer :: ps => null()
      type(hdf_data), pointer :: psd => null()
      integer , parameter :: N_PARAMS = 4
      integer , parameter :: N_PARAMS_EMS = 2
      integer :: i , last , first
      
     
      if ( .not. file_test ( self % file )) then 
         print*, 'file not available channel' 
         stop
      end if
   
      allocate ( sds_name ( N_PARAMS) )
      sds_name = (/ character (len =20) :: 'albedo' , 'transmission' , 'spherical_albedo', 'reflectance'  /)
      
      if ( hdf_get_file_sds ( self%file, nsds , sds , nsdsn = N_PARAMS, sds_name = sds_name ) < 0 ) then
         print*,'hdf file not readable'
         stop
      end if   
      deallocate ( sds_name )
      
      
      call self % alloc 
      ps => sds(1); psd=> ps%data
      self % cld_alb = reshape(psd%r4values,[9,29,45])
      
      ps => sds(2); psd=>ps%data
      self % cld_trn =reshape (psd%r4values,[9,29,45])
      
      ps => sds(3); psd=> ps%data
      self % cld_sph_alb = reshape ( psd%r4values, [9,29] )
      
      ps => sds(4); psd=>ps%data
      
      ! reconstruct 5d array
      do i = 1 , 45 
         first = (i -1 ) *  9 * 29 * 45 * 45 + 1
         last = first + (9 * 29 * 45 * 45) - 1
         self % cld_refl(:,:,:,:,i) = reshape( psd%r4values(first : last  ),[9,29,45,45])
      end do 
       
      if ( self % has_ems ) then
         allocate ( sds_name_ems ( N_PARAMS_EMS) )
         sds_name_ems = (/ character(len=20) :: 'cloud_emissivity' , 'cloud_transmission' /)
      
         if ( hdf_get_file_sds ( self%file_ems, nsds , sds , nsdsn = N_PARAMS_EMS, sds_name = sds_name_ems ) < 0 ) then
            print*,'hdf file not readable'
            stop
         end if   
         deallocate ( sds_name_ems )
         ps => sds(1); psd=> ps%data
         self % cld_ems = reshape(psd%r4values,[9,29,45])
         
         ps => sds(2); psd=> ps%data
         self % cld_trn_ems = reshape(psd%r4values,[9,29,45])
         
         
      end if
      
      ps => null()
      psd => null()
      
     
   end subroutine lut_data__read_hdf
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__alloc ( self)
      class ( lut_data_type ) :: self
      
      allocate ( self % cld_alb (9,29,45))   
      allocate ( self % cld_trn (9,29,45))
      allocate ( self % cld_sph_alb (9,29))
      allocate ( self % cld_refl (9,29,45,45,45))
      allocate ( self % cld_ems (9,29,45))
      allocate ( self % cld_trn_ems (9,29,45))
   end subroutine lut_data__alloc
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__dealloc ( self)
      class ( lut_data_type ) :: self
      
      if ( allocated (self % cld_alb) ) deallocate ( self % cld_alb )   
      if ( allocated (self % cld_trn) ) deallocate ( self % cld_trn )
      if ( allocated (self % cld_sph_alb) ) deallocate ( self % cld_sph_alb )
      if ( allocated (self % cld_refl) ) deallocate ( self % cld_refl )
      if ( allocated (self % cld_ems) ) deallocate ( self % cld_ems )
      if ( allocated (self % cld_trn_ems) ) deallocate ( self % cld_trn_ems )
      
      self % is_set = .false.
      
   end subroutine lut_data__dealloc   
   
   
    
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__init_dims( self)        
      class ( lut_type ) :: self
      character ( len = 300 )  :: hdf_file
      integer :: nsds                     

      character(len=MAXNCNAM), dimension(:), allocatable :: &
       sds_name

      type(hdf_sds), dimension(:), allocatable, target :: &
       sds                          ! Tableau des structures des SDS extraits

      integer :: isds
      integer  :: i                                      
       
      type(hdf_sds), pointer                           :: &
       ps                           ! Pointeur sur la structure du SDS courant 
      
      type(hdf_data), pointer                          :: &
       psd                        ! Pointeur sur les données du SDS courant
            
      
           
      hdf_file = self % channel ( 1) % phase (1 ) % file
    
		if ( .not. file_test(hdf_file) ) then
			print*,'lut file not existing! ==> ', trim(hdf_file)
			stop
		end if
		
      ! this should be read from file, but this is also possible
      self %  dims% n_sat_zen = 45
      self %  dims% n_sol_zen = 45
      self %  dims% n_rel_azi = 45
      self %  dims% n_cod = 29
      self %  dims% n_cps = 9
      
      call  self % dims % dealloc
      call  self % dims % alloc 
      
      allocate ( sds_name ( 5)) 
      sds_name =(/ character(len=30) :: 'sensor_zenith_angle'  &
            , 'solar_zenith_angle' &
            , 'relative_azimuth_angle' &
            , 'log10_optical_depth' &
            , 'log10_eff_radius'/)
      
       
      if (hdf_get_file_sds(hdf_file, nsds, sds, nsdsn = 5, sds_name = sds_name) < 0) then
         print*,'hdf file not readable ', trim(hdf_file)
         stop  
      end if   
      ps => sds(1); psd=> ps%data
      self %  dims% sat_zen = psd%r4values 
       
      ps => sds(2); psd=> ps%data
      self %  dims% sol_zen = psd%r4values 
      
      ps => sds(3); psd=> ps%data
      self %  dims% rel_azi = psd%r4values 
      
      ps => sds(4); psd=> ps%data
      self %  dims% cod = psd%r4values 
      
      ps => sds(5); psd=> ps%data      
      self %  dims% cps = psd%r4values 
      
      psd => null()
      
      deallocate ( sds_name)
            
   end subroutine lut__init_dims
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_dim__alloc ( self )
      class ( lut_dim_type) :: self
      
      allocate ( self % sat_zen (self % n_sat_zen)  )
      allocate ( self % sol_zen (self % n_sol_zen)  )
      allocate ( self % rel_azi (self % n_rel_azi)  )
      allocate ( self % cod (self % n_cod)  )
      allocate ( self % cps (self % n_cps)  )
      
      
   
   end subroutine lut_dim__alloc
  
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_dim__dealloc ( self )
      class ( lut_dim_type) :: self
      
      if ( allocated (self % sat_zen) ) deallocate ( self % sat_zen  )
      if ( allocated (self % sol_zen) ) deallocate ( self % sol_zen   )
      if ( allocated (self % rel_azi) ) deallocate ( self % rel_azi  )
      if ( allocated (self % cod ) ) deallocate ( self % cod )
      if ( allocated (self % cps ) ) deallocate ( self % cps  )
   end subroutine lut_dim__dealloc
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__getProperty ( self, sat , sol , azi , sensor )
      class ( lut_type) :: self
      real , intent(out) , optional :: sat
      real , intent(out) , optional :: sol   
      real , intent(out) , optional :: azi
      character(10), intent(out), optional :: sensor
      
      if ( present(sensor)) sensor = self % sensor
      if ( present(sat) ) sat = self % dims % sat_zen ( self % pos_sat )
      if ( present(sol) ) sol = self % dims % sat_zen ( self % pos_sol )
      if ( present(azi) ) azi = self % dims % sat_zen ( self % pos_azi )
   
   end subroutine lut__getProperty
   

end module  dcomp_lut_mod
