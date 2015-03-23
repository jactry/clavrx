! $Id$
module dcomp_interface_TYPEs_mod

! --  Module works as interface between CLAVR-x and DCOMP
!
!

   ! - parameters
   integer, parameter :: REAL4 = selected_real_kind(6,37)
   integer, parameter :: INT1 = selected_int_kind(1)   
   integer, parameter :: INT2 = selected_int_kind(2)
   real ( kind = real4), parameter :: MISSING_REAL4 = -999.
   
   ! - CLAVRX uses 44 MODIS/VIIRS channels
   integer , parameter :: N_CHN = 44   
   
   ! - object for 2D real4 arrays
   type d2_real4_type
      logical :: is_set
      integer :: xdim
      integer :: ydim
      real ( kind = real4 ) , dimension(:,:) , allocatable  :: d  
     
   end type d2_real4_type
   
   ! - object for 2D int1 arrays
   type d2_int1_type
      logical :: is_set
      integer :: xdim
      integer :: ydim
      integer ( kind = int1 ) , dimension(:,:) , allocatable  :: d  
     
   end type d2_int1_type
   
   ! - object for 2D int2 arrays
   type d2_int2_type
      logical :: is_set
      integer :: xdim
      integer :: ydim
      integer ( kind = int2 ) , dimension(:,:) , allocatable  :: d  
   end type d2_int2_type
   
   ! - object for 2D logical arrays
   type d2_flag_type
      logical :: is_set
      logical, dimension(:,:), allocatable :: d
     
   end type d2_flag_type
   
   ! - object for gas coeff values
   type gas_coeff_type
      logical :: is_set
      real ( kind = real4 )  :: d (3)  
   end type  gas_coeff_type
   
   
   ! - main dcomp input type
   type dcomp_in_type
   
      ! - configure
      integer :: mode
      character ( len = 1024) :: lut_path
      integer :: sensor_wmo_id
      logical :: is_channel_on (N_CHN)
      
      ! - satellite input
      TYPE ( d2_real4_TYPE) ,allocatable :: refl(:)!(N_CHN)
      TYPE ( d2_real4_TYPE) , allocatable :: rad (:)!(N_CHN)
      TYPE ( d2_real4_TYPE) :: sat
      TYPE ( d2_real4_TYPE) :: sol
      TYPE ( d2_real4_TYPE) :: azi
      
      ! - cloud products
      TYPE ( d2_int1_TYPE )  :: cloud_mask
      TYPE ( d2_real4_TYPE ) :: cloud_type
      TYPE ( d2_real4_TYPE ) :: cloud_hgt
      TYPE ( d2_real4_TYPE ) :: cloud_temp
      TYPE ( d2_real4_TYPE ) :: cloud_press
      TYPE ( d2_real4_TYPE ) :: tau_acha
      
      ! - flags
      TYPE ( d2_flag_TYPE )  :: is_land 
      TYPE ( d2_flag_TYPE )  :: is_valid
      
      
      ! - surface
      TYPE ( d2_real4_TYPE) ,allocatable :: alb_sfc (:)!( N_CHN )
      TYPE ( d2_real4_TYPE ) :: press_sfc
      TYPE ( d2_real4_TYPE) ,allocatable :: emiss_sfc (:)!( N_CHN )
      TYPE ( d2_int1_TYPE)   :: snow_class 
      
      ! - atmosphere
      TYPE ( d2_real4_TYPE ) :: ozone_nwp
      TYPE ( d2_real4_TYPE ) :: tpw_ac
      TYPE ( d2_real4_TYPE ) ,allocatable:: trans_ac_nadir (:)!( N_CHN )
      TYPE ( d2_real4_TYPE ) ,allocatable:: rad_clear_sky_toc(:)! ( N_CHN )
      TYPE ( d2_real4_TYPE ),allocatable :: rad_clear_sky_toa(:)! ( N_CHN )
      
      ! - coeffecients,params
      real :: sun_earth_dist
      TYPE ( gas_coeff_type )  :: gas_coeff ( N_CHN)
      real :: solar_irradiance(N_CHN)
      
      contains
      final :: in_destructor 
         
   end type dcomp_in_type
   
   
   ! - DCOMP output
      
   type dcomp_out_type
      type ( d2_real4_type) :: cod  
      type ( d2_real4_type) :: cps
      type ( d2_real4_type) :: cod_unc
      type ( d2_real4_type) :: ref_unc
      type ( d2_real4_type) :: cld_trn_sol
      type ( d2_real4_type) :: cld_trn_obs
      type ( d2_real4_type) :: cld_alb
      type ( d2_real4_type) :: cld_sph_alb
      type ( d2_int2_type) :: quality
      type ( d2_int2_type) :: info
      type ( d2_real4_type) :: iwp  
      type ( d2_real4_type) :: lwp

	   character ( len = 200 ) :: version
	   real :: successrate
      integer :: nr_clouds
      integer :: nr_obs
      integer :: nr_success_cod
      integer :: nr_success_cps
      
      contains
      final :: out_destructor
   end type dcomp_out_type 
   
   ! - Enumerated cloud type
   type et_cloud_type_type
      integer(kind=int1) :: FIRST = 0
      integer(kind=int1) :: CLEAR = 0
      integer(kind=int1) :: PROB_CLEAR = 1
      integer(kind=int1) :: FOG = 2
      integer(kind=int1) :: WATER = 3
      integer(kind=int1) :: SUPERCOOLED = 4
      integer(kind=int1) :: MIXED = 5
      integer(kind=int1) :: OPAQUE_ICE = 6
      integer(kind=int1) :: TICE = 6
      integer(kind=int1) :: CIRRUS = 7
      integer(kind=int1) :: OVERLAP = 8
      integer(kind=int1) :: OVERSHOOTING = 9
      integer(kind=int1) :: UNKNOWN = 10
      integer(kind=int1) :: DUST = 11
      integer(kind=int1) :: SMOKE = 12
      integer(kind=int1) :: FIRE = 13  
      integer(kind=int1) :: LAST = 13
   end type
   
   type ( et_cloud_type_type ) , protected :: EM_cloud_type
   
   ! - Enumerated clod mask
   type et_cloud_mask_type
      integer(kind=int1) :: LAST = 3
      integer(kind=int1) :: CLOUDY = 3
      integer(kind=int1) :: PROB_CLOUDY = 2
      integer(kind=int1) :: PROB_CLEAR = 1
      integer(kind=int1) :: CLEAR = 0
      integer(kind=int1) :: FIRST = 0
   end type
   
   type ( et_cloud_mask_type ) , protected :: EM_cloud_mask
   
   ! - Enumerated snow/sea ice class
   type  et_snow_class_type
      integer(kind=int1) :: FIRST = 1
      integer(kind=int1) :: NO_SNOW = 1
      integer(kind=int1) :: SEA_ICE = 2
      integer(kind=int1) :: SNOW = 3
      integer(kind=int1) :: LAST = 3
   end type 
   
   type (  et_snow_class_type ) , protected :: EM_snow_class 
   
   interface alloc_dcomp
      module procedure &
      alloc_it_d2_real, alloc_it_d2_int, alloc_it_d2_log
   
   end interface 
   
   interface dcomp_out_type
      module procedure new_output
   end interface dcomp_out_type
   
   interface dcomp_in_type
      module procedure new_input
   end interface dcomp_in_type
   
   !interface alloc_dcomp
   !   subroutine alloc_it_d2_real ( str , xdim , ydim )
   !      import   d2_real4_type
   !      type ( d2_real4_type ) :: str
   !      integer :: xdim , ydim     
   !   end subroutine alloc_it_d2_real
      
   !   subroutine alloc_it_d2_int ( str , xdim , ydim )
   !      import d2_int1_type
   !      type ( d2_int1_type ) :: str
   !      integer :: xdim , ydim 
   !   end subroutine alloc_it_d2_int
      
   !   subroutine alloc_it_d2_log ( str , xdim , ydim )
   !      import d2_flag_type
   !      type ( d2_flag_type ) :: str
   !      integer :: xdim , ydim 
   !   end subroutine alloc_it_d2_log
   !end interface alloc_dcomp 
   
contains
   
   
   function new_output ( dim_1, dim_2 )
      integer , intent(in) :: dim_1, dim_2
      type ( dcomp_out_type ) :: new_output
      
      allocate ( new_output % cod % d         ( dim_1 , dim_2))
      allocate ( new_output % cps % d         ( dim_1 , dim_2))
      allocate ( new_output % cod_unc % d     ( dim_1 , dim_2))
      allocate ( new_output % ref_unc % d     ( dim_1 , dim_2))
      allocate ( new_output % cld_trn_sol % d ( dim_1 , dim_2))
      allocate ( new_output % cld_trn_obs % d ( dim_1 , dim_2))
      allocate ( new_output % cld_alb % d     ( dim_1 , dim_2))
      allocate ( new_output % cld_sph_alb % d ( dim_1 , dim_2))
      allocate ( new_output % info % d        ( dim_1 , dim_2))
      allocate ( new_output % quality % d     ( dim_1 , dim_2))
      allocate ( new_output % lwp % d         ( dim_1 , dim_2))
      allocate ( new_output % iwp % d         ( dim_1 , dim_2))
      
      new_output % cod % d           =  MISSING_REAL4 
      new_output % cps % d           =  MISSING_REAL4 
      new_output % cod_unc % d       =  MISSING_REAL4
      new_output % ref_unc % d       =  MISSING_REAL4 
      new_output % cld_trn_sol % d   =  MISSING_REAL4  
      new_output % cld_trn_obs % d   =  MISSING_REAL4   
      new_output % cld_alb % d       =  MISSING_REAL4  
      new_output % cld_sph_alb % d   =  MISSING_REAL4 
  
   end function new_output
   
   
   function new_input ( dim_1, dim_2, chan_on )
      integer, intent(in) :: dim_1, dim_2
      logical, intent(in) :: chan_on ( :)
      type ( dcomp_in_type ) :: new_input
      integer :: n_chn
      integer :: idx_chn
      
         ! === ALLOCATION
         
      allocate (  new_input % refl (N_CHN) ) 
      allocate (  new_input % rad (N_CHN) )
      
      allocate ( new_input %  alb_sfc (N_CHN) )
      allocate ( new_input %  emiss_sfc (N_CHN) )
      allocate ( new_input %  rad_clear_sky_toa (N_CHN) ) 
      allocate ( new_input %  rad_clear_sky_toc (N_CHN) ) 
      allocate ( new_input %  trans_ac_nadir (N_CHN) ) 
      
         
      n_chn = size ( chan_on)   
      do idx_chn = 1 , n_chn
        
         if ( .not. chan_on (idx_chn) ) cycle
         
         new_input % is_channel_on(idx_chn) = .true.
         call  alloc_dcomp ( new_input % refl    (  idx_chn  ) , dim_1,dim_2 ) 
         call  alloc_dcomp ( new_input % alb_sfc (  idx_chn  ) ,  dim_1,dim_2 ) 
                 
         if ( idx_chn >= 20 ) then   
            call  alloc_dcomp ( new_input % rad (  idx_chn  ) ,  dim_1,dim_2 )  
            call  alloc_dcomp ( new_input % emiss_sfc (  idx_chn  ) ,  dim_1,dim_2 )
            call  alloc_dcomp ( new_input % rad_clear_sky_toa ( idx_chn ),  dim_1,dim_2 )
            call  alloc_dcomp ( new_input % rad_clear_sky_toc ( idx_chn ), dim_1,dim_2 )
            call alloc_dcomp (  new_input % trans_ac_nadir ( idx_chn )   , dim_1,dim_2 )
         end if   
      end do
      
      call  alloc_dcomp ( new_input % snow_class,   dim_1,dim_2 )
      call  alloc_dcomp ( new_input % is_land,      dim_1,dim_2 )
      call  alloc_dcomp ( new_input % cloud_press,  dim_1,dim_2 )
      call  alloc_dcomp ( new_input % cloud_temp,   dim_1,dim_2 )
      call  alloc_dcomp ( new_input % cloud_hgt,    dim_1,dim_2 )
      call  alloc_dcomp ( new_input % cloud_type,   dim_1,dim_2 )
      call  alloc_dcomp ( new_input % cloud_mask,   dim_1,dim_2 )
      call  alloc_dcomp ( new_input % tau_acha,   dim_1,dim_2 )
      call  alloc_dcomp ( new_input % ozone_nwp,    dim_1,dim_2 )
      call  alloc_dcomp ( new_input % tpw_ac,       dim_1,dim_2 )
      call  alloc_dcomp ( new_input % press_sfc,    dim_1,dim_2 )
      call  alloc_dcomp ( new_input % is_valid,     dim_1,dim_2 )
      call  alloc_dcomp ( new_input % sol,          dim_1,dim_2 )
      call  alloc_dcomp ( new_input % sat,          dim_1,dim_2 )
      call  alloc_dcomp ( new_input % azi,          dim_1,dim_2 )
      
     
   
   
   
   end function new_input
   
   subroutine out_destructor ( this )
      type ( dcomp_out_type) :: this
      
	   if (allocated ( this % cod % d ) ) deallocate (  this % cod % d ) 
		if ( allocated ( this % cps % d )) deallocate ( this % cps % d)
	   if (allocated ( this % cod_unc % d ) ) deallocate (  this % cod_unc % d ) 
		if ( allocated ( this % ref_unc % d )) deallocate ( this % ref_unc % d)   
      if (allocated ( this % cld_trn_sol % d ) ) deallocate (  this % cld_trn_sol % d ) 
		if ( allocated ( this % cld_trn_obs % d )) deallocate ( this % cld_trn_obs % d)
	   if (allocated ( this % cld_alb % d ) ) deallocate (  this % cld_alb % d ) 
		if ( allocated ( this % cld_sph_alb % d )) deallocate ( this % cld_sph_alb % d)
		if (allocated ( this % quality % d ) ) deallocate (  this % quality % d ) 
		if ( allocated ( this % info % d )) deallocate ( this % info % d)
	   if (allocated ( this % iwp % d ) ) deallocate (  this % iwp % d ) 
		if ( allocated ( this % lwp % d )) deallocate ( this % lwp % d)
   end subroutine out_destructor
   
   
   !  --  allocation routines   
   subroutine alloc_it_d2_real ( str , xdim , ydim )
      type ( d2_real4_type ) :: str
      integer :: xdim , ydim
      integer :: alloc_stat = 0
         
      allocate ( str % d ( xdim,  ydim) , stat = alloc_stat)
      if ( alloc_stat /= 0 ) then
         print*,'alloc error'            
      end if      
   end subroutine alloc_it_d2_real
   
   !  --  allocation routine
   subroutine alloc_it_d2_int ( str , xdim , ydim )
      type ( d2_int1_type ) :: str
      integer :: xdim , ydim
      integer :: alloc_stat = 0
         
      allocate ( str % d ( xdim,  ydim) , stat = alloc_stat)
      if ( alloc_stat /= 0 ) then
         print*,'alloc error'            
      end if      
   end subroutine alloc_it_d2_int
   
   !  --  allocation routine   
   subroutine alloc_it_d2_log ( str , xdim , ydim )
      type ( d2_flag_type ) :: str
      integer :: xdim , ydim
      integer :: alloc_stat = 0
         
      allocate ( str % d ( xdim,  ydim) , stat = alloc_stat)
      if ( alloc_stat /= 0 ) then
         print*,'alloc error'            
      end if      
   end subroutine alloc_it_d2_log
   
   
   
   
   
   subroutine in_destructor ( this )
      type ( dcomp_in_type ) :: this
      integer :: i
        
        
      do i = 1, N_CHN 
         if ( allocated (this % refl(i) % d) ) deallocate ( this % refl(i) % d )
			if ( allocated (this % alb_sfc(i) % d) ) deallocate ( this % alb_sfc(i) % d )
         if ( allocated (this % rad(i) % d) ) deallocate ( this % rad(i) % d )
         if ( allocated (this % emiss_sfc(i) % d) ) deallocate ( this % emiss_sfc(i) % d )
         if ( allocated (this % trans_ac_nadir(i) % d) ) deallocate ( this % trans_ac_nadir(i) % d )
			if ( allocated (this % rad_clear_sky_toa(i) % d) ) deallocate ( this % rad_clear_sky_toa(i) % d )
			if ( allocated (this % rad_clear_sky_toc(i) % d) ) deallocate ( this % rad_clear_sky_toc(i) % d ) 
                    
      end do
      
      if ( allocated (this % refl) ) deallocate ( this % refl )
      if ( allocated (this % alb_sfc) ) deallocate ( this % alb_sfc )
      if ( allocated (this % rad) ) deallocate ( this % rad )
      if ( allocated (this % emiss_sfc) ) deallocate ( this % emiss_sfc )
      if ( allocated (this % trans_ac_nadir) ) deallocate ( this % trans_ac_nadir )
      if ( allocated (this % rad_clear_sky_toa) ) deallocate ( this % rad_clear_sky_toa )
      if ( allocated (this % rad_clear_sky_toc) ) deallocate ( this % rad_clear_sky_toc ) 
      
      
      if ( allocated (this % sol % d) ) deallocate ( this % sol  % d )
      if ( allocated (this % sat % d) ) deallocate ( this % sat  % d )
      if ( allocated (this % azi % d) ) deallocate ( this % azi  % d )
      
      if ( allocated (this % cloud_mask % d) )  deallocate ( this % cloud_mask  % d )
      if ( allocated (this % cloud_type % d) )  deallocate ( this % cloud_type  % d )
      if ( allocated (this % cloud_hgt % d) )   deallocate ( this % cloud_hgt % d )
      if ( allocated (this % cloud_temp % d) )  deallocate ( this % cloud_temp  % d )
      if ( allocated (this % cloud_press % d) ) deallocate ( this % cloud_press  % d )
      if ( allocated (this % tau_acha % d) )    deallocate ( this % tau_acha  % d )
      
      if ( allocated (this % snow_class % d) ) deallocate ( this % snow_class  % d )
      if ( allocated (this % is_land % d) ) deallocate ( this % is_land  % d )
    	if ( allocated (this % ozone_nwp % d) )    deallocate ( this %  ozone_nwp   % d )
		if ( allocated (this % tpw_ac % d) )    deallocate ( this %  tpw_ac   % d )
		if ( allocated (this % press_sfc % d) )    deallocate ( this %  press_sfc % d )
		if ( allocated (this % is_valid % d) ) deallocate ( this % is_valid  % d )
      
      
      
   end subroutine in_destructor
   
   
   
      
   
end module dcomp_interface_types_mod
