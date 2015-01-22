!   $Id:$
!
!
!   
module cx_rtm_mod

   use cx_nwp_mod, only: nwp_main_type
  
   use cx_geo_mod, only: &
      geo_type
      
   use cx_rtm_tools_mod,only: &
      nlevels_rtm &
        , emissivity &
        , convert_profiles_nwp_rtm 
         
   use cx_tools_science_planck_mod, only: &
      planck_bt2rad &
      , planck_rad2bt
   
   use sfc_data_mod, only: &
      sfc_main_type
      
  
      
   use cx_tools_science_mod, only : &
       locate
   use cx_sat_mod, only: &
      sat_main_type
   
   
   
   real, parameter :: Pi = 3.14159
   real, parameter :: DTOR = PI / 180.
   real, parameter, private:: CO2_RATIO = 380.0 !in ppmv
   
   
    integer, parameter :: real4 = selected_real_kind(6,37)
   integer, parameter :: int1 = selected_int_kind(1)
   
   integer, parameter, public:: RTM_NVZEN = 50
   real, parameter, public::  RTM_VZA_BINSIZE = 1./ RTM_NVZEN
   
  
   !---------------------------------------------------------------------
   ! RTM structure definition
   !---------------------------------------------------------------------
   
   !  +++++++   NWP Grid   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! ---   NWP Grid Profile
   TYPE, public :: rtm_ngrid_chn_prof_type
      logical :: is_set = .false.
      REAL (kind=real4), DIMENSION(:), allocatable :: Trans_Two_Way_Atm_Clr
      REAL (kind=real4), DIMENSION(:), allocatable :: Trans_Prof
      REAL (kind=real4), DIMENSION(:), allocatable :: Trans_Atm_Clr_Solar_total
      REAL (kind=real4), DIMENSION(:), allocatable :: Trans_Atm_Clr_solar
      REAL (kind=real4), DIMENSION(:), allocatable :: Rad_Prof
      REAL (kind=real4), DIMENSION(:), allocatable :: Rad_BB_Prof
      contains
      procedure :: allocate => allocate_rtm_ngrid_chn_prof
   end TYPE Rtm_ngrid_chn_prof_type
    
   ! ---  NWP Grid channel 
   type :: rtm_ngrid_chn_type
      logical :: is_set = .false.
      type (rtm_ngrid_chn_prof_type ),  allocatable :: chn (:)
      contains
      procedure :: allocate => allocate_rtm_ngrid_chn
   end type rtm_ngrid_chn_type
  
   ! ---  NWP ngrid sub
   type, public :: rtm_ngrid_sub_type
      logical :: is_set = .false.
      integer  :: Level_sfc = 0
      integer  :: Level_tropo = 0
      integer  :: Inversion_Level = 0
      integer  :: Level440 = 0
      integer  :: Level850 = 0
      TYPE (Rtm_ngrid_chn_type), allocatable :: angl (:)
      
      REAL (kind=real4), dimension(:), allocatable :: T_Prof
      REAL (kind=real4), dimension(:), allocatable :: Z_Prof
      REAL (kind=real4), dimension(:), allocatable :: Wvmr_Prof
      REAL (kind=real4), dimension(:), allocatable :: Ozmr_Prof
      contains
      procedure :: allocate => allocate_rtm_ngrid_sub
   end type rtm_ngrid_sub_type
   
   !  +++++++   Satellite Grid   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! --- sgrid sub channel data structure  ------------
   type rtm_sgrid_chn_data_type
      logical :: is_set
      REAL (kind=real4) :: rad_atm        ! -- radiance at toa emitted by atmosphere
      REAL (kind=real4) :: trans_atm      ! -- transmission of atmosphere from surface to toa 
      REAL (kind=real4) :: rad_atm_sfc    ! -- radiance at toa emitted by atmosphere and surface
      REAL (kind=real4) :: bt_atm_sfc     ! -- bt computes with planck from rad_atm_sfc
      real ( kind = real4) :: emiss_tropo ! -- emissivity at tropopause
      real ( kind = real4 ) :: emiss_clear  ! - emiss at toa , (or?)
      real ( kind = real4 ) :: refl_clear ! - assumed visible reflectance ( channels 1-20)
   end type rtm_sgrid_chn_data_type
   
   ! --- sgrid sub channel  structure   -------
   type rtm_sgrid_chn_type
      logical :: is_set
      type (rtm_sgrid_chn_data_type) , allocatable :: chn(:)
      integer :: idx_angl
      integer :: level_sfc
      contains
      procedure :: allocate => allocate_sgrid_chn
   end type rtm_sgrid_chn_type
   
   !  +++++++   MAIN   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! --- RTM main  --------------- 
   type rtm_main_type
      logical :: is_set  
      logical :: is_sgrid_set 
      type(rtm_ngrid_sub_type), public, allocatable :: ngrid (:,:)
      type(rtm_sgrid_chn_type), public, allocatable :: sgrid (:,:)
      contains
      procedure :: allocate_ngrid
      procedure :: deallocate_ngrid
      procedure :: allocate_sgrid
      procedure :: deallocate_sgrid
      procedure :: populate
      procedure :: set_key_levels
      
   end type rtm_main_type
       
   integer , public:: pfaast_chn_idx (36)
  
 
contains
   ! ---
   !  called in process_clavrx.f90
   !  --
   subroutine populate ( self,  nwp , geo , sfc,  sat  )
   
      implicit none
      class (rtm_main_type) , target :: self
      type(nwp_main_type ) :: nwp
      type(geo_type) , target :: geo
      type(sfc_main_type) :: sfc
      real , pointer, dimension(:,:)  :: satzen
      
      type ( sat_main_type ) , target :: sat
      
      
      integer :: i_rtm
      integer :: j_rtm
      integer :: i_sat
      integer :: j_sat
      
      
      integer :: min_idx_x 
      integer :: max_idx_x 
      integer :: min_idx_y 
      integer :: max_idx_y 
     
      real, allocatable :: rad_sfc_bb(:,:,:)
      
      integer :: num_sat_pxls
      integer :: Idx_Ang
      real :: Satzen_Mid_bin
      character ( len = 128) :: Ancil_Data_Dir
     
      
      integer :: kbin
      integer :: chan_idx
      
      real :: t_midlev 
      real :: t_lev
      real :: rad_midlev 
      real :: rad_lev
      real, pointer :: trans_prof(:)
      real, pointer :: t_prof_pnt(:)
      real, pointer :: rad_prof(:)
      real, pointer :: rad_bb_prof(:)
      integer :: k
      integer :: idx_ang_min,idx_ang_max
      real :: sfc_rad
      
      integer, pointer :: i_nwp_x , i_nwp_y, i_angl 
      integer :: i_sfc , i_tropo
      
      type (rtm_ngrid_chn_prof_type) , pointer :: ngrid_p
      type (rtm_sgrid_chn_data_type) , pointer :: sgrid_p 
      
      real :: rad_ch20_temp
      real :: sfc_ref_ch20
      real :: rad_refl_factor
      real :: get_rad_refl_factor
      real :: solar_radiance_ch20
      
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      satzen => sat % geo % sat_zen
      
      ! - create mapping of pfaast channel numbers to sensor channels
      call map_pfaast_chn
      
      ! -  loop over all rtm       
      min_idx_x = minval ( geo % idx_nwp_x , geo % idx_nwp_x > 0)
      max_idx_x = maxval ( geo % idx_nwp_x , geo % idx_nwp_x > 0)
      min_idx_y = minval ( geo % idx_nwp_y , geo % idx_nwp_y > 0)
      max_idx_y = maxval ( geo % idx_nwp_y , geo % idx_nwp_y > 0)
      
      print*,min_idx_x,max_idx_x,min_idx_y,max_idx_y
      
     
      call self % allocate_ngrid( nwp % n_lon, nwp % n_lat)
      print*, nwp % n_lon, nwp % n_lat
      
         
      rtm_x_loop: do i_rtm = min_idx_x, max_idx_x
         num_sat_pxls = count (  geo % idx_nwp_x == i_rtm )
         if ( num_sat_pxls == 0) cycle rtm_x_loop
         
         rtm_y_loop: do j_rtm =  min_idx_y, max_idx_y
            num_sat_pxls = count ( geo % idx_nwp_y == j_rtm .and. geo % idx_nwp_x == i_rtm )
           
            if ( num_sat_pxls == 0) cycle rtm_y_loop
              
               !  allocate the atmospheric profiles
               call self % ngrid (i_rtm , j_rtm ) % allocate( nlevels_rtm)
                             
               ! - transform to RTM profiles
               ! t_nwp, z_nwp, ozmr_nwp, rh_nwp
               !TODO make abetter input and output strutucre for self routine!
              
               call convert_profiles_NWP_RTM  ( &
                     &   nwp % p_std & 
                     & , nwp % sfc_level(i_rtm, j_rtm) &
                     & , nwp % psfc (i_rtm , j_rtm ) &
                     & , nwp % wvmrsfc (i_rtm, j_rtm) &
                     & , nwp % tmpair (i_rtm, j_rtm) &
                     & , nwp % t_prof (i_rtm, j_rtm, :) &  
                     & , nwp % ozone_prof (i_rtm, j_rtm, :) &
                     & , nwp % wvmr_prof (i_rtm, j_rtm, :) &
                     & , nwp % z_prof (i_rtm, j_rtm,:) &
                      
                     &  , self % ngrid (i_rtm, j_rtm) % z_prof &
                     &  , self % ngrid (i_rtm, j_rtm) % t_prof &
                     &  , self % ngrid (i_rtm, j_rtm) % wvmr_prof &
                     &  , self % ngrid (i_rtm, j_rtm) % ozmr_prof &
                     )
                     
                     
             
               ! allocate and fill profiles in rtm structure
                             
               idx_ang_min = floor(50*minval(cos ( satzen * DTOR) &
                              , geo % idx_nwp_y == j_rtm .and. &
                               geo % idx_nwp_x == i_rtm ))
               
               idx_ang_max =  ceiling(50* maxval(cos (satzen * DTOR) &
                               , geo % idx_nwp_y == j_rtm .and. &
                               geo % idx_nwp_x == i_rtm))
               
              
               rtm_ang_loop: do Idx_ang = idx_ang_min , idx_ang_max
                  
                  ! - check if self was already computed ( is self possible?)
                  if ( self % ngrid (i_rtm, j_rtm) % angl(idx_ang) % is_set) cycle rtm_ang_loop
                  
                  self % ngrid (i_rtm, j_rtm) % angl(idx_ang) % is_set = .true.
                  
                  ! 
                  call self % ngrid (i_rtm, j_rtm) % angl(idx_ang) % allocate()
                  
                  
                  !--- determine the zenith angle for the RTM
                  Satzen_Mid_Bin = acos (  DTOR *((Idx_Ang-1)*Rtm_vza_binsize + Rtm_vza_binsize/2.0))  
                  Ancil_Data_Dir = '/DATA/Ancil_Data/clavrx_ancil_data/'
                
               
                  rtm_chn_loop: do Chan_Idx = 20 , 36
                     if (  sat % chan_on (chan_idx) .eqv. .false. ) cycle
                     kbin = pfaast_chn_idx(chan_Idx)
                     if (kbin == 0 ) cycle rtm_chn_loop
                     
                     call self % ngrid (i_rtm, j_rtm) % angl (Idx_ang) % chn (chan_idx) % allocate( nlevels_rtm)  
                     
                     
                    
                     t_prof_pnt => self % ngrid (i_rtm, j_rtm) % t_prof
                     trans_prof => self % ngrid (i_rtm, j_rtm) % angl (Idx_ang) % chn (chan_idx) % trans_prof
                     rad_prof => self % ngrid (i_rtm, j_rtm) % angl (Idx_ang) % chn (chan_idx) % rad_prof
                     rad_bb_prof => self % ngrid (i_rtm, j_rtm) % angl (Idx_ang) % chn (chan_idx) % rad_bb_prof
                     
                     ! pfaast
                     !TODO for all sensors ...
                     
                     call tran_viirsm(Ancil_Data_Dir &
                        , self % ngrid (i_rtm, j_rtm) % T_Prof &
                        , self % ngrid (i_rtm, j_rtm) % Wvmr_Prof &
                        , self % ngrid (i_rtm, j_rtm) % Ozmr_Prof &
                        , Satzen_Mid_Bin &
                        , CO2_RATIO &
                        , kbin &
                        , trans_prof )
                        
                                          
                     level_loop: do k = 2 , nLevels_rtm
                        
                     
                     !    acronym lev stands for representative value for self level ( the mean between k and k+1)
                                          
                        t_midlev = 0.5 * (  t_prof_pnt(k-1) + t_prof_pnt (k) )  
                       
                        call planck_bt2rad ( t_midlev , 'VIIRS' , chan_idx , rad_midlev )
                        
                        ! - see planck comment
                                              
                        rad_prof(k) = rad_prof(k-1) + (trans_prof(k-1) - trans_prof(k)) * rad_midlev
                        
                        t_lev = t_prof_pnt (k)
                        call planck_bt2rad ( t_lev , 'VIIRS' , chan_idx , rad_lev )
                        
                         rad_bb_prof(k) = rad_prof(k) +  trans_prof(k) * rad_lev
                         
                                                                               
                     end do level_loop
                     
                        
                     trans_prof => null()   
                     t_prof_pnt => null()
                     rad_prof => null()
                     rad_bb_prof => null()
                    
                     
                  end do rtm_chn_loop  
               end do   rtm_ang_loop  
            end do rtm_y_loop
         end do   rtm_x_loop
     
      
      ! set key levels
      call self % set_key_levels (  nwp )
       
      !
      !  - pixel stuff
      !  - pixels has specific
      !     surface idx
      
      !
      ! 1. allocate
     
      call self % allocate_sgrid( geo % n_x, geo % n_y)
      
      ! - compute black body BT from surface
      allocate ( rad_sfc_bb (20:36, geo % n_x, geo % n_y ))
      
      do chan_idx = 20, 36
         call planck_bt2rad ( nwp % sgrid % t_sfc , 'VIIRS' , chan_idx, rad_sfc_bb(chan_idx,:,:) )
      end do
      
      self % sgrid  % idx_angl = floor( 50 * cos (satzen * DTOR ))
     
      ! - loop over sat grid
      do i_sat = 1 ,  geo % n_x
         do j_sat = 1,  geo % n_y
           
            i_nwp_x => geo%idx_nwp_x(i_sat,j_sat)
            i_nwp_y => geo%idx_nwp_y(i_sat,j_sat)
            i_angl => self % sgrid(i_sat,j_sat)  % idx_angl
            if ( i_nwp_x < 0 .or. i_nwp_y < 0 ) cycle
            ! - determine sfc level in profile
            ! initial set to the coarse nwp grid sfc
            i_sfc = self % ngrid (i_nwp_x, i_nwp_y) % level_sfc 
           
            call locate ( self % ngrid (i_nwp_x, i_nwp_y) % z_prof ,nlevels_rtm &
               , (real(sfc % elevation % data (i_sat,j_sat) ) )/1000.0 , i_sfc )
            
            i_tropo = self % ngrid (i_nwp_x, i_nwp_y) % level_tropo
            
            call self % sgrid (i_sat,j_sat)  % allocate()
            do chan_idx = 20, 36
               
               ngrid_p => self % ngrid (i_nwp_x, i_nwp_y) % angl (I_angl)  % chn(Chan_idx)
                  
               sgrid_p => self % sgrid(i_sat,j_sat)  % chn(Chan_idx)
           
               if (  sat % chan_on (chan_idx) .eqv. .false. ) cycle
             
               if ( .not. ngrid_p % is_set ) cycle
              
                                      
               sgrid_p % rad_atm = ngrid_p % Rad_Prof (i_sfc)
               sgrid_p % trans_atm =  ngrid_p %   Trans_Prof (i_sfc)
               !-TODO weights for exact 
                        
               sfc_rad = sfc % emis ( chan_idx) % data (i_sat,j_sat)  * rad_sfc_bb(chan_idx,i_sat , j_sat )  
              
               sgrid_p % rad_atm_sfc = &
                   sgrid_p % rad_atm + sgrid_p % trans_atm * sfc_rad
              
               ! - add solar part to self mixed signal channel    
               if ( chan_idx == 20 ) then
                  
                  sfc_ref_ch20 = ( 1 - sfc % emis ( chan_idx) % data (i_sat,j_sat) ) 
                  rad_refl_factor = get_rad_refl_factor ( 'VIIRS' , sat % geo % sol_zen (i_sat , j_sat) )
                  
                  solar_radiance_ch20 = sgrid_p % trans_atm &
                     & * sfc_ref_ch20  &
                     & *  1 / rad_refl_factor
                  
                  sgrid_p % rad_atm_sfc = &
                     & sgrid_p % rad_atm_sfc &
                     & + solar_radiance_ch20
                       
               end if    
                 
                       
               call planck_rad2bt ( sgrid_p % rad_atm_sfc &
                         , 'VIIRS' &
                         , chan_idx &
                         , sgrid_p % bt_atm_sfc )  
                        
               self % ngrid (i_nwp_x, i_nwp_y) % angl (I_angl) % chn (chan_idx) % is_set = .true.   
               
               
               if ( sat % chn (chan_idx) % rad (i_sat,j_sat) > 0. .and. sfc % emis ( chan_idx) % data (i_sat,j_sat) > 0 ) then
                  
                  sgrid_p % emiss_tropo = emissivity ( sat % chn (chan_idx) % rad (i_sat,j_sat) &
                                          & , sgrid_p % rad_atm_sfc &
                                          & , ngrid_p % Rad_bb_Prof (i_tropo))
                 
                                    
               end if
            end do
           
         
            
            ! - channel 20 radiance from bt108 
            call planck_bt2rad ( self % sgrid(i_sat,j_sat)  % chn(31) % bt_atm_sfc , 'VIIRS' , 20, rad_ch20_temp )
           
            self % sgrid(i_sat,j_sat)  % chn(20)  % emiss_clear = &
               & self % sgrid(i_sat,j_sat)  % chn(20) %rad_atm_sfc / rad_ch20_temp  
           
            ! - solar clear sky for channel 1
            !TODO atmospheric corrrection
            self % sgrid(i_sat , j_sat) % chn(1) % refl_clear = sfc % modis_w_sky(1) % data (i_sat,j_sat)  
             
         end do
      end do
 
      i_nwp_x => null()
      i_nwp_y => null()
      i_angl => null()
      
      sgrid_p => null()
      ngrid_p => null()
     
      
      deallocate ( rad_sfc_bb)
 
       !call rtm_solar (sfc , rtm )
      
   end subroutine populate

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
   !  ++++++++++ RTM Main structure  ++++++++++++++
   !   -------  NGRID ALLOC
   subroutine allocate_ngrid ( self , n_lon , n_lat )
      class ( rtm_main_type ) , target :: self
      integer :: n_lon , n_lat
      integer :: astatus
            
      if ( .not. self % is_set ) then 
         allocate ( self % ngrid (n_lon , n_lat), stat = astatus)
         if ( astatus /= 0 ) then
            print "(a,'Not enough memory to allocate Rtm_Params structure.')"
         end if
         
         
         self % is_set = .true.
         
        
      end if   
   end subroutine allocate_ngrid

   !  -------    NGRID DEALLOC   ------
   subroutine deallocate_ngrid (self)
      class ( rtm_main_type ) :: self
      
      if (self % is_set ) then
         if (allocated (self % ngrid ) ) deallocate ( self % ngrid ) 
      end if
      self % is_set = .false.   
   end subroutine deallocate_ngrid
  
   !  -----------  SGRID ALLOC   ------
   subroutine allocate_sgrid ( self , n_x , n_y )
      class ( rtm_main_type ) :: self
      integer :: n_x , n_y
      integer :: astatus
            
      if ( .not. self % is_sgrid_set ) then 
         allocate ( self % sgrid (n_x , n_y), stat = astatus)
         if ( astatus /= 0 ) then
            print "(a,'Not enough memory to allocate satellite grid Rtm_Params structure.')"
         end if
       
         self % is_sgrid_set = .true.
      end if   
   end subroutine allocate_sgrid
   
   ! --------------  SGRID DEALLOC  ------
   subroutine deallocate_sgrid (self)
      class ( rtm_main_type ) :: self
     
      if (self % is_sgrid_set ) then
         if (allocated (self % sgrid)) deallocate ( self % sgrid) 
      end if
      self % is_sgrid_set = .false.
         
   end subroutine deallocate_sgrid
   
   ! +++++++ NGRID Substrctures ++++++
   subroutine allocate_rtm_ngrid_sub ( self , n_levels)
      class (rtm_ngrid_sub_type) :: self
      integer :: n_levels
     
      
      if ( .not. self%is_set ) then 
         allocate ( self % t_prof(n_levels))
         allocate ( self % z_prof(n_levels))
         allocate ( self % wvmr_prof(n_levels))
         allocate ( self % ozmr_prof(n_levels))
         allocate ( self % angl (RTM_NVZEN) )
      end if
      self%is_set = .true.
   end subroutine allocate_rtm_ngrid_sub
   
   !  -----------------------------------------------------
   subroutine allocate_rtm_ngrid_chn (self)
      class (rtm_ngrid_chn_type) :: self 
      allocate ( self % chn(36))   
   end subroutine allocate_rtm_ngrid_chn
     
   !  -----------------------------------------------------
   subroutine allocate_rtm_ngrid_chn_prof ( self , n_levels)
      class (rtm_ngrid_chn_prof_type) :: self
      integer :: n_levels
      
      allocate ( self % Trans_Prof ( n_levels ) )
      allocate ( self % Rad_Prof ( n_levels ) )
      allocate ( self % Rad_BB_Prof ( n_levels ) )
      allocate ( self % Trans_Two_Way_Atm_Clr( n_levels ) )
      self % is_set = .true.
         
   end subroutine allocate_rtm_ngrid_chn_prof
   
   ! +++++++ SGRID Substrctures ++++++
   !  -----------------------------------------------------
   subroutine allocate_sgrid_chn (self)
      class (rtm_sgrid_chn_type) :: self 
      if (.not. allocated (self % chn) ) allocate ( self % chn(36))   
   end subroutine allocate_sgrid_chn
  
  
    !
   !
   !
   subroutine set_key_levels (this, nwp)
      
      use cx_nwp_mod
      use cx_rtm_constants_mod, only : &
          nlevels_rtm &
          , p_std_rtm
      
      class ( rtm_main_type) :: this
      type ( nwp_main_type ) :: nwp
      
      integer :: k_nwp , i_nwp, j_nwp
      
      ! -- ngrid levels
      this % ngrid  % level_sfc = nlevels_rtm
      
      do i_nwp = 1 , nwp % n_lon 
         
         do j_nwp = 1, nwp % n_lat 
        
            do k_nwp = NLevels_Rtm , 1 , -1
               if ( p_std_rtm ( k_nwp) < nwp   % psfc ( i_nwp, j_nwp))then
                  this % ngrid ( i_nwp, j_nwp) % level_sfc = k_nwp
                  exit
               end if
            end do
            
            do k_nwp = 1 , this % ngrid ( i_nwp, j_nwp) % level_sfc - 1
            
               if ( p_std_rtm (k_nwp) <= nwp  % p_trop( i_nwp, j_nwp) .and. &
                   p_std_rtm (k_nwp+1) > nwp  % p_trop( i_nwp, j_nwp) ) then
                   this % ngrid ( i_nwp, j_nwp) % level_tropo = k_nwp
               end if  
                         
            end do 
            
           
            
         end do
      end do   
      
      
      
      !--------------------------------------------------------------------
      !--- find tropopause Level  based on tropopause pressure
      !--- tropopause is between tropopause_Level and tropopaue_Level + 1
      !--------------------------------------------------------------------
      !do k_nwp = 1, Rtm_grid % Sfc_Level-1
      !   if ((P_Std_Rtm(k_nwp) <= nwp % P_trop (Lon_Idx,Lat_Idx) ) .and. &
      !      (P_Std_Rtm(k+1) > P_trop_Nwp(Lon_Idx,Lat_Idx))) then
      !      Rtm_grid % Tropo_Level = k
      !   end if
      !end do
      
      
      
      ! -- sgrid 
      
   
   end subroutine set_key_levels
   
   subroutine rtm_solar ( sfc , rtm ) 
      use sfc_data_mod, only: sfc_main_type
      use cx_rtm_types_mod
      type(sfc_main_type) :: sfc
      type (rtm_main_type) , target :: rtm
      
      !-TODO atmosph correction
      rtm % sgrid(22,22) % chn(1) % refl_clear = sfc % modis_w_sky(1) % data (22,22)
      
   end subroutine rtm_solar
   
   !
   !
   
   subroutine map_pfaast_chn
       ! for viirs  
      pfaast_chn_idx = 0
      pfaast_chn_idx (20) = 1
      pfaast_chn_idx (22) = 2
      pfaast_chn_idx (29) = 3
      pfaast_chn_idx (31) = 4
      pfaast_chn_idx (32) = 5
   
   end subroutine map_pfaast_chn

  

end module cx_rtm_mod
