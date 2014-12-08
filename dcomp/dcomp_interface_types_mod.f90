! $Id$
module dcomp_interface_TYPEs_mod

! --  Module works as interface between CLAVR-x and DCOMP
!
!

   ! - parameters
   integer, parameter :: REAL4 = selected_real_kind(6,37)
   integer, parameter :: INT1 = selected_int_kind(1)   
   integer, parameter :: INT2 = selected_int_kind(2)
   
   ! - CLAVRX uses 42 MODIS/VIIRS channels
   integer , parameter :: N_CHN = 42   
   
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
      real ( kind = real4 ) , dimension(3) :: d   
   end type  gas_coeff_type
   
   
   ! - main dcomp input type
   type dcomp_in_type
   
      ! - configure
      integer :: mode
      character ( len = 1024) :: lut_path
      integer :: sensor_wmo_id
      logical :: is_channel_on (N_CHN)
      
      ! - satellite input
      TYPE ( d2_real4_TYPE) :: refl(N_CHN)
      TYPE ( d2_real4_TYPE) :: rad (N_CHN)
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
      TYPE ( d2_real4_TYPE)  :: alb_sfc (N_CHN)
      TYPE ( d2_real4_TYPE ) :: press_sfc
      TYPE ( d2_real4_TYPE)  :: emiss_sfc(N_CHN)
      TYPE ( d2_int1_TYPE)   :: snow_class 
      
      ! - atmosphere
      TYPE ( d2_real4_TYPE ) :: ozone_nwp
      TYPE ( d2_real4_TYPE ) :: tpw_ac
      TYPE ( d2_real4_TYPE ) :: trans_ac_nadir ( N_CHN )
      TYPE ( d2_real4_TYPE ) :: rad_clear_sky_toc ( N_CHN )
      TYPE ( d2_real4_TYPE ) :: rad_clear_sky_toa ( N_CHN )
      
      ! - coeffecients,params
      real :: sun_earth_dist
      TYPE ( gas_coeff_type ), dimension(40) :: gas_coeff
      real :: solar_irradiance(40)
         
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
   
   
   subroutine deallocate_dcompin ( dcomp_str )
      type ( dcomp_in_type ) , intent (inout) :: dcomp_str
      integer :: i
        
        
      do i = 1, N_CHN 
         if ( allocated (dcomp_str % refl(i) % d) ) deallocate ( dcomp_str % refl(i) % d )
			if ( allocated (dcomp_str % alb_sfc(i) % d) ) deallocate ( dcomp_str % alb_sfc(i) % d )
         if ( allocated (dcomp_str % rad(i) % d) ) deallocate ( dcomp_str % rad(i) % d )
         if ( allocated (dcomp_str % emiss_sfc(i) % d) ) deallocate ( dcomp_str % emiss_sfc(i) % d )
         if ( allocated (dcomp_str % trans_ac_nadir(i) % d) ) deallocate ( dcomp_str % trans_ac_nadir(i) % d )
			if ( allocated (dcomp_str % rad_clear_sky_toa(i) % d) ) deallocate ( dcomp_str % rad_clear_sky_toa(i) % d )
			if ( allocated (dcomp_str % rad_clear_sky_toc(i) % d) ) deallocate ( dcomp_str % rad_clear_sky_toc(i) % d )
			
            
      end do
      if ( allocated (dcomp_str % sol % d) ) deallocate ( dcomp_str % sol  % d )
      if ( allocated (dcomp_str % sat % d) ) deallocate ( dcomp_str % sat  % d )
      if ( allocated (dcomp_str % azi % d) ) deallocate ( dcomp_str % azi  % d )
      
      if ( allocated (dcomp_str % cloud_mask % d) )  deallocate ( dcomp_str % cloud_mask  % d )
      if ( allocated (dcomp_str % cloud_type % d) )  deallocate ( dcomp_str % cloud_type  % d )
      if ( allocated (dcomp_str % cloud_hgt % d) )   deallocate ( dcomp_str % cloud_hgt % d )
      if ( allocated (dcomp_str % cloud_temp % d) )  deallocate ( dcomp_str % cloud_temp  % d )
      if ( allocated (dcomp_str % cloud_press % d) ) deallocate ( dcomp_str % cloud_press  % d )
      if ( allocated (dcomp_str % tau_acha % d) )    deallocate ( dcomp_str % tau_acha  % d )
      
      if ( allocated (dcomp_str % snow_class % d) ) deallocate ( dcomp_str % snow_class  % d )
      if ( allocated (dcomp_str % is_land % d) ) deallocate ( dcomp_str % is_land  % d )
    	if ( allocated (dcomp_str % ozone_nwp % d) )    deallocate ( dcomp_str %  ozone_nwp   % d )
		if ( allocated (dcomp_str % tpw_ac % d) )    deallocate ( dcomp_str %  tpw_ac   % d )
		if ( allocated (dcomp_str % press_sfc % d) )    deallocate ( dcomp_str %  press_sfc % d )
		if ( allocated (dcomp_str % is_valid % d) ) deallocate ( dcomp_str % is_valid  % d )
   end subroutine deallocate_dcompin
   
   
   subroutine deallocate_dcompout (dcomp_str_out)
		type ( dcomp_out_type) , intent (inout) :: dcomp_str_out
		
		
	   if (allocated ( dcomp_str_out % cod % d ) ) deallocate (  dcomp_str_out % cod % d ) 
		if ( allocated ( dcomp_str_out % cps % d )) deallocate ( dcomp_str_out % cps % d)
	   if (allocated ( dcomp_str_out % cod_unc % d ) ) deallocate (  dcomp_str_out % cod_unc % d ) 
		if ( allocated ( dcomp_str_out % ref_unc % d )) deallocate ( dcomp_str_out % ref_unc % d)   
      if (allocated ( dcomp_str_out % cld_trn_sol % d ) ) deallocate (  dcomp_str_out % cld_trn_sol % d ) 
		if ( allocated ( dcomp_str_out % cld_trn_obs % d )) deallocate ( dcomp_str_out % cld_trn_obs % d)
	   if (allocated ( dcomp_str_out % cld_alb % d ) ) deallocate (  dcomp_str_out % cld_alb % d ) 
		if ( allocated ( dcomp_str_out % cld_sph_alb % d )) deallocate ( dcomp_str_out % cld_sph_alb % d)
		if (allocated ( dcomp_str_out % quality % d ) ) deallocate (  dcomp_str_out % quality % d ) 
		if ( allocated ( dcomp_str_out % info % d )) deallocate ( dcomp_str_out % info % d)
	   if (allocated ( dcomp_str_out % iwp % d ) ) deallocate (  dcomp_str_out % iwp % d ) 
		if ( allocated ( dcomp_str_out % lwp % d )) deallocate ( dcomp_str_out % lwp % d)
      
 
	end subroutine deallocate_dcompout
   
      
   
end module dcomp_interface_types_mod
