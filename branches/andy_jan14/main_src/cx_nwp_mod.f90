! $Header:$
! 
!   01/21/2015: in clavrx
!
module cx_nwp_mod
   use date_tools_mod
    
   use cr_physics_mod, only: &
     wvmr_from_rh_temp_press
     
   implicit none
   
    integer, parameter :: real4 = selected_real_kind(6,37)
   
   type nwp_sgrid_type
      logical :: is_set
      real ( kind = real4 ) :: t_sfc
      real ( kind = real4 ) :: t_tropopause
      real ( kind = real4 ) :: t_air
      real ( kind = real4 ) :: rh_sfc
      real ( kind = real4 ) :: p_sfc
      real ( kind = real4 ) :: p_msl ! what is msl?
      real ( kind = real4 ) :: tpw
      real ( kind = real4 ) :: ozone
      real ( kind = real4 ) :: k_index
     
      real ( kind = real4 ) :: sc_lwp  !  what is this?
      real ( kind = real4 ) :: lwp
      real ( kind = real4 ) :: iwp
      real ( kind = real4 ) :: cwp
      real ( kind = real4 ) :: pc
      real ( kind = real4 ) :: cloud_fraction
      real ( kind = real4 ) :: ncld_layers
      real ( kind = real4 ) :: cld_type
      real ( kind = real4 ) :: wind_speed
      real ( kind = real4 ) :: wind_dir
      real ( kind = real4 ) :: lifting_condensation_lev
      real ( kind = real4 ) :: convection_condensation_lev
      
        
   end type  nwp_sgrid_type
   
   
   type , public :: nwp_main_type
      integer, dimension(:,:), allocatable, public :: Mask
      integer, dimension(:,:), allocatable, public :: bad_mask
      integer, dimension(:,:), allocatable, public :: sfc_type
      real , dimension(:,:), allocatable, public :: Satzen
      real , dimension(:,:), allocatable, public :: Solzen
      real , dimension(:,:), allocatable, public :: Psfc
      real , dimension(:,:), allocatable, public :: Pmsl
      real , dimension(:,:), allocatable, public :: Zsfc
      real , dimension(:,:), allocatable, public :: Tmpsfc
      real , dimension(:,:), allocatable, public :: Tmpsfc_Before
      real , dimension(:,:), allocatable, public :: Tmpsfc_After
      real , dimension(:,:), allocatable, public :: Tmpair
      real , dimension(:,:), allocatable, public :: Tmpair_uni
      real , dimension(:,:), allocatable, public :: Tmpsfc_uni
      real , dimension(:,:), allocatable, public :: T_Trop
      real , dimension(:,:), allocatable, public :: P_Trop
      real , dimension(:,:), allocatable, public :: Rhsfc
      real , dimension(:,:), allocatable, public :: wvmrsfc
      real , dimension(:,:), allocatable, public :: Tpw
      real , dimension(:,:), allocatable, public :: Uth
      real , dimension(:,:), allocatable, public :: Hght500
      real , dimension(:,:), allocatable, public :: Ozone
      real , dimension(:,:), allocatable, public :: weasd
      real , dimension(:,:), allocatable, public :: U_Wnd_10m
      real , dimension(:,:), allocatable, public :: V_Wnd_10m
      real , dimension(:,:), allocatable, public :: Wnd_Spd_10m
      real , dimension(:,:), allocatable, public :: Wnd_Dir_10m
      integer , dimension(:,:), allocatable, public :: Land
      real , dimension(:,:), allocatable, public :: Ice
      integer , dimension(:,:), allocatable, public :: Tropo_Level
      integer , dimension(:,:), allocatable, public :: Level850
      integer , dimension(:,:), allocatable, public :: Level700
      integer , dimension(:,:), allocatable, public :: Level500
      integer , dimension(:,:), allocatable, public :: Sfc_Level
      integer , dimension(:,:), allocatable, public :: Inversion_Level
      integer , dimension(:,:,:), allocatable, public :: Inversion_Level_Profile
      integer , dimension(:,:), allocatable, public :: Strato_Level
      real , dimension(:,:), allocatable, public :: Lifting_Condensation_Level_Height !km
      real , dimension(:,:), allocatable, public :: Convective_Condensation_Level_Height !km
      real , dimension(:,:), allocatable, public :: Freezing_Level_Height !km
      real , dimension(:,:), allocatable, public :: Upper_Limit_Water_Height !km
      real , dimension(:,:), allocatable, public :: K_Index !K
      real , dimension(:,:), allocatable, public :: Pc
      real , dimension(:,:), allocatable, public :: Lwp
      real , dimension(:,:), allocatable, public :: Iwp
      real , dimension(:,:), allocatable, public :: Cwp
      real , dimension(:,:), allocatable, public :: Cloud_Fraction_Satellite
      real , dimension(:,:), allocatable, public :: High_Cloud_Fraction_Satellite
      real , dimension(:,:), allocatable, public :: Mid_Cloud_Fraction_Satellite
      real , dimension(:,:), allocatable, public :: Low_Cloud_Fraction_Satellite
      integer , dimension(:,:), allocatable, public :: Ncld_Layers
      integer , dimension(:,:), allocatable, public :: Cld_Type
      real , dimension(:,:), allocatable, public :: Lat
      real, dimension(:,:), allocatable, public :: Lon
      real, public :: lat_first
      real, public :: lon_first
      integer , public :: d_lat
      integer , public :: d_lon
      integer , public :: n_lat
      integer , public :: n_lon      
      real  , dimension(:), allocatable, public :: P_Std
      real , dimension(:,:,:), allocatable, public :: Z_Prof
      real  , dimension(:,:,:), allocatable, public :: T_Prof
      real , dimension(:,:,:), allocatable, public :: Rh_Prof
      real  , dimension(:,:,:), allocatable, public :: Ozone_Prof
	   real  , dimension(:,:,:), allocatable, public :: wvmr_prof
      real , dimension(:,:,:), allocatable, public :: Tpw_Prof
      real , dimension(:,:,:), allocatable, public :: Clwmr_Prof
      real , dimension(:,:,:), allocatable, public :: U_Wnd_Prof
      real , dimension(:,:,:), allocatable, public :: V_Wnd_Prof
      real , dimension(:,:), allocatable, public :: temp2d_1
      real , dimension(:,:), allocatable, public :: temp2d_2
      real , dimension(:,:,:), allocatable, public :: temp3d_1
      real , dimension(:,:,:), allocatable, public :: temp3d_2
      real , dimension(:,:,:), allocatable, public :: temp3d
  integer , public :: n_levels
   
      integer, public:: start_hour
      integer, public:: end_hour
      
      
      type ( nwp_sgrid_type) , allocatable :: sgrid (:,:) 
   
   contains
      procedure :: populate => populate_all_nwp
      procedure :: deallocate =>deallocate_all
      procedure :: set_special_levels
      procedure :: assign_to_sat_grid 
      procedure :: deallocate_sat_grid
   
   end type nwp_main_type

contains
   ! ----------------------------------------------------------------
   ! 
   ! -----------------------------------------------------------------
   subroutine populate_all_nwp ( this , date_start , date_end , ancil_path , nwp_opt)
      
           
      implicit none
      class ( nwp_main_type) :: this
      type ( date_type ) , intent (in) :: date_start , date_end
      character (256) , intent(in) ::ancil_path
      integer , intent(in) :: nwp_opt
   
      call read_gfs_data ( this, date_start , date_end , ancil_path, nwp_opt ) 
      call set_special_levels( this  )
          
   end subroutine populate_all_nwp
   
   
   !  -----------------------------------------
   !   set_special_levels
   !   Andi Walther Nov 2013
   !  ----------------------------------------------------------  ------------------------
   !      find special key levels in NWP grid and profile levels
   !      2013/12/30 AW
   !      should be in a science sub tool, becuase this is not a pure read routine
   !  ----------------------------------------------------------  --------------------------
   subroutine set_special_levels ( this )
      implicit none
      class ( nwp_main_type) :: this
      integer :: i_nwp, j_nwp, k_nwp
      
      integer :: level850, level700, level500
      real :: p_trop_temp
      
      allocate ( this % sfc_level (this % n_lon, this % n_lat) ) 
      allocate ( this % tropo_level (this % n_lon, this % n_lat) ) 
      
      this % sfc_level = 0
      
      
      do i_nwp = 1 , this % n_lon 
         do j_nwp = 1 , this % n_lat
            
            ! - surface level
            do k_nwp = this % n_levels , 1 , -1 
               if ( this % p_std (k_nwp) < this % psfc ( i_nwp,j_nwp) ) then
                  this % sfc_level (i_nwp,j_nwp) = k_nwp
                  
                  exit
               end if   
            end do
            
            ! - standard levels needed only here ( probably)
            do k_nwp = 1 , this % sfc_level (i_nwp,j_nwp) -1
               if ( ( this % p_std ( k_nwp) <= 850.0 ) .and. ( this % p_std ( k_nwp) > 850.0 )) then
                  level850 = k_nwp 
               end if
               if ( ( this % p_std ( k_nwp) <= 700.0 ) .and. ( this % p_std ( k_nwp) > 700.0 )) then
                  level700 = k_nwp 
               end if
               if ( ( this % p_std ( k_nwp) <= 500.0 ) .and. ( this % p_std ( k_nwp) > 500.0 )) then
                  level500 = k_nwp 
               end if
            
            end do
            
            ! - tropopause  level  based on tropopause pressure
            !--- tropopause is between tropopause_level and tropopaue_level + 1
            !--------------------------------------------------------------------
            
            !--- constrain tropopause pressure to be greater than 75 mb
            p_trop_temp = max ( this % p_trop ( i_nwp,j_nwp) , 75.0)
            this % tropo_level ( i_nwp,j_nwp)  = 1
            do k_nwp = 1 , this % sfc_level (i_nwp,j_nwp) -1
               if ( ( this % p_std (k_nwp) <= p_trop_temp )  &
                  .and. ( this % p_std (k_nwp + 1) > p_trop_temp ) ) then
                  this % tropo_level ( i_nwp,j_nwp)  = k_nwp
               end if
            end do
            
         end do
      end do   
      

   end subroutine set_special_levels
   
   !  -----------------------------------------
   !   READ_GFS_DATA
   !   Andi Walther Nov 2013
   !  -------------------------------------------------
   subroutine read_gfs_data(gNWP , date_start , date_end , ancil_path , nwp_opt)
      
      use file_tools, only: &
         file_test &
      , file_basename
      
      ! use hdf
      !  use cr_physics_mod
      implicit none
      
      type ( date_type ) , intent (in) :: date_start , date_end
      type ( nwp_main_type), intent (in out)  :: gNWP
      character (256) , intent(in) ::ancil_path
      integer, intent(in) :: nwp_opt
      
      character(len = 255 ) :: nwp_file1, nwp_file2
      type ( date_type ) :: after_gfs_time
      type ( date_type ) :: before_gfs_time

      ! - file and hdf data handling
      integer :: sd_id_1, sd_id_2
      integer :: istatus
      integer :: sfstart, sfend, sfrnatt, sffattr, sfrcatt, sfrdata, sfendacc

      ! - nwp attributes
      integer :: gfs_n_levels
      integer :: gfs_n_levels_o3
      integer :: gfs_n_levels_rh
      integer :: gfs_n_levels_clw
      integer :: gfs_n_lat
      integer :: gfs_n_lon
	   real :: gfs_d_lat
	   real :: gfs_d_lon
	   real :: gfs_lat_1
	   real :: gfs_lon_1
      character (3) :: gfs_array_order_1 ,  gfs_array_order_2

      integer, dimension ( 2 ) :: sds_start_2d, sds_stride_2d, sds_edge_2d
      integer :: sfselect , sds_id , sfn2index
	   real :: time_wgt
      integer ::  ii ,  jj , i , j

	   real , dimension(:), allocatable, target :: P_Std_Nwp
      integer :: DFACC_READ =1 
      
      ! - gfs are 12 h forecasts !
      before_gfs_time =  next_6h ( date_start , - 3 )
      after_gfs_time =  next_6h ( date_start , -2 )
      
      ! -  the minus 2 because this is a 12h forecast to stay between 0 and 1.
      time_wgt = time_diff_weight (date_start, before_gfs_time, after_gfs_time) -2.
     
       select case (nwp_opt)
       case (1)  
         nwp_file1 = trim(trim(ancil_path)//'dynamic/gfs'//"/"// &
                      "gfs." // before_gfs_time % yymmddhh //"_F012.hdf")
         nwp_file2 = trim(trim(ancil_path)//'dynamic/gfs'//"/"// &
                      "gfs." // after_gfs_time % yymmddhh //"_F012.hdf")
                      
        case (2)
            print*,'baeeh ncep later'
        
        case (3)
              nwp_file1 = trim(trim(ancil_path)//'dynamic/cfsr'//"/"// &
                      "gfs." // before_gfs_time % yymmddhh //"_F012.hdf")
            nwp_file2 = trim(trim(ancil_path)//'dynamic/cfsr'//"/"// &
                      "gfs." // after_gfs_time % yymmddhh //"_F012.hdf")
        
        case(4)
              nwp_file1 = trim(trim(ancil_path)//'../gdas'//"/"// &
                      "gfs." // before_gfs_time % yymmddhh //"_F012.hdf")
            nwp_file2 = trim(trim(ancil_path)//'../gdas'//"/"// &
                      "gfs." // after_gfs_time % yymmddhh //"_F012.hdf")
        
        end select              
     
      nwp_file1 = trim(nwp_file1)
      if ( .not. file_test (nwp_file1)) then 
    
         print*,'file missing nwp1 => ',nwp_file1
         stop
      end if
  
      nwp_file2 = trim(nwp_file2)

      if ( .not. file_test (nwp_file2)) then 
    
         print*,'file missing nwp2 ..',nwp_file2
         stop
      end if
   
      ! - check if files exist and open					
      if ((file_test (nwp_file1) .eqv. .false. ) .or. (file_test (nwp_file1) .eqv. .false. )) then
         print*, 'nwp files are not existing ', file_basename (nwp_file1) 
      end if
      

      sd_id_1 = sfstart ( nwp_file1, DFACC_READ)
      sd_id_2 = sfstart ( nwp_file2, DFACC_READ)

      ! - reads attribute data from file 1
      istatus = 0
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF PRESSURE LEVELS"),gfs_n_levels) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF O3MR LEVELS"),gfs_n_levels_o3) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF RH LEVELS"),gfs_n_levels_rh) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF CLWMR LEVELS"),gfs_n_levels_clw) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF LATITUDES"),gfs_n_lat) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF LONGITUDES"),gfs_n_lon) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"LATITUDE RESOLUTION"),gfs_d_lat) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"LONGITUDE RESOLUTION"),gfs_d_lon) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"FIRST LATITUDE"),gfs_lat_1) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"FIRST LONGITUDE"),gfs_lon_1) + istatus
      istatus = sfrcatt(sd_id_1,sffattr(sd_id_1,"3D ARRAY ORDER"), gfs_array_order_1) + istatus

   ! - don't read file2 attributes ( assume they are the same) 
  istatus = sfrcatt(sd_id_2,sffattr(sd_id_2,"3D ARRAY ORDER"), gfs_array_order_2) + istatus
     
      gNWP % d_lat  = gfs_d_lat
   gNWP % d_lon = gfs_d_lon
      
      
    
      
      !-TODO
   
  
      allocate (gNWP%psfc(gfs_n_lon, gfs_n_lat))
   allocate (gNWP%lon(gfs_n_lon, gfs_n_lat))
   allocate (gNWP%lat(gfs_n_lon, gfs_n_lat))

   allocate (gNWP%pmsl(gfs_n_lon, gfs_n_lat))
   allocate (gNWP%tmpsfc(gfs_n_lon, gfs_n_lat))
   allocate (gNWP%tmpair(gfs_n_lon, gfs_n_lat))
      allocate (gNWP%zsfc(gfs_n_lon, gfs_n_lat))
   allocate (gNWP%rhsfc(gfs_n_lon, gfs_n_lat))
   allocate (gNWP%land(gfs_n_lon, gfs_n_lat))
   allocate (gNWP%weasd(gfs_n_lon, gfs_n_lat))
      allocate (gNWP%ice(gfs_n_lon, gfs_n_lat))
      allocate (gNWP%t_trop(gfs_n_lon, gfs_n_lat))
      allocate (gNWP%p_trop(gfs_n_lon, gfs_n_lat))
      allocate (gNWP%ozone(gfs_n_lon, gfs_n_lat))
      allocate (gNWP%tpw(gfs_n_lon, gfs_n_lat))
    
      do i = 1 , gfs_n_lon 
         do j = 1 , gfs_n_lat 
            gNWP % lon ( i , j ) =  gfs_lon_1 + ( i - 1 ) * gfs_d_lon
            gNWP % lat ( i , j ) =  gfs_lat_1  + ( j -1 ) * gfs_d_lat
         end do
      end do
      
      gNWP % lat_first = gfs_lat_1
      gNWP % lon_first = gfs_lon_1
     
      gNWP % ice    = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'ice fraction', time_wgt )
      gNWP % psfc   = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'surface pressure', time_wgt )
      gNWP % land   = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'land mask', time_wgt )

      gNWP % pmsl   = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'MSL pressure', time_wgt )
      gNWP % tmpsfc = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'surface temperature', time_wgt )
      
      gNWP % zsfc   = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'surface height', time_wgt )
      gNWP % tmpair = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'temperature at sigma=0.995', time_wgt )
      gNWP % rhsfc  = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'rh at sigma=0.995', time_wgt )
     
      gNWP % weasd  = read_nwp_2d  ( sd_id_1 , sd_id_2 , 'water equivalent snow depth', time_wgt )
      gNWP % t_trop = read_nwp_2d ( sd_id_1, sd_id_2 , 'tropopause temperature', time_wgt ) 
      gNWP % p_trop = read_nwp_2d ( sd_id_1, sd_id_2 , 'tropopause pressure', time_wgt )
      gNWP % ozone = read_nwp_2d ( sd_id_1, sd_id_2 , 'total ozone', time_wgt )
      gNWP % tpw = read_nwp_2d ( sd_id_1, sd_id_2 , 'total precipitable water', time_wgt )

      istatus =0
      allocate ( gNWP%p_std ( gfs_n_levels) )
      sds_id = sfselect(sd_id_1, sfn2index(sd_id_1,'pressure levels'))
      
         
      Istatus = sfrdata(sds_id,  0 , 1 , gfs_n_levels  ,  &
                      gNWP%p_std) + Istatus
      Istatus = sfendacc(sds_id) + Istatus 
     
      gNWP%n_levels = gfs_n_levels

      allocate ( gNWP%t_prof(gfs_n_lon, gfs_n_lat,gfs_n_levels))
      allocate ( gNWP%z_prof(gfs_n_lon, gfs_n_lat,gfs_n_levels))
      allocate ( gNWP%rh_prof(gfs_n_lon, gfs_n_lat,gfs_n_levels))
      allocate ( gNWP%ozone_prof(gfs_n_lon, gfs_n_lat,gfs_n_levels))
      allocate ( gNWP%wvmr_prof(gfs_n_lon, gfs_n_lat,gfs_n_levels))

      gNWP%wvmr_prof = 0.

      gNWP%t_prof = read_nwp_3d (sd_id_1 , sd_id_2 , 'temperature' , time_wgt )
      gNWP%z_prof = read_nwp_3d (sd_id_1 , sd_id_2 , 'height' , time_wgt )
      gNWP%rh_prof = read_nwp_3d (sd_id_1 , sd_id_2 , 'rh' , time_wgt )
      gNWP%ozone_prof = read_nwp_3d (sd_id_1 , sd_id_2 , 'o3mr' , time_wgt )

      do ii =1 , gfs_n_lon
         do jj =1 , gfs_n_lat
            gNWP%wvmr_prof(ii,jj,:) =  wvmr_from_rh_temp_press (  gNWP % rh_prof ( ii , jj , : )  &
                , gNWP % t_prof ( ii , jj , : ) , gNWP % p_std(:) )
         end do
      end do 
 
      allocate (gNWP%wvmrsfc(gfs_n_lon, gfs_n_lat))
      !gNWP % wvmrsfc =  wvmr_from_rh_temp_press (  gNWP % rhsfc  &
      !                , gNWP % tmpair , gNWP % pmsl )

      deallocate (gNWP%rh_prof)
      istatus = sfend(sd_id_1)
      istatus = sfend(sd_id_2)
      
      
      gNWP % n_lat = gfs_n_lat
      gNWP % n_lon = gfs_n_lon
      
      ! -  fix gfs RH bug
  
      ! - other scientific computations ( wind speed , tropopasue etc)
    
      ! - be happy
    
   contains
   
      ! ======================================================
      !
      ! ======================================================
      
   function read_nwp_2d (sd_id_1 , sd_id_2, sds_name, time_wgt)  &
                      &   result (result_real)
      implicit none
 
         integer, intent(in):: sd_id_1, sd_id_2
      character(len=*), intent(in):: sds_name
	      real, intent(in) :: time_wgt
      integer, dimension(2) :: sds_start, sds_stride, sds_edge
	      real , dimension(:,:), allocatable :: result_real
  
      integer :: sds_id_1, sds_id_2
	      real , dimension(:,:), allocatable :: dum1, dum2
      integer :: istatus_1 = 0
         integer :: istatus_2 = 0
         integer :: istatus = 0
      integer :: gfs_n_lon , gfs_n_lat
 
      ! - hdf functions
      integer :: sfselect ,  sfn2index, sdreaddata
      integer :: sffattr
      
         !  executable   
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF LATITUDES"),gfs_n_lat) + istatus
         istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF LONGITUDES"),gfs_n_lon) + istatus
  
      sds_start = [ 0, 0 ]

      sds_stride = [ 1, 1 ]
      sds_edge = [ gfs_n_lon , gfs_n_lat ]

      Istatus_1 = 0

         sds_id_1 = sfselect(sd_id_1, sfn2index( sd_id_1, sds_name ) )
      sds_id_2 = sfselect(sd_id_2, sfn2index( sd_id_2, sds_name ) )
  
      !--- read sds's
      allocate( dum1  ( gfs_n_lon , gfs_n_lat ) , dum2( gfs_n_lon , gfs_n_lat ))
      allocate(result_real (gfs_n_lon , gfs_n_lat ))
  
         Istatus_1 = sfrData(sds_id_1,sds_start, sds_stride, sds_edge, dum1) + Istatus_1
         Istatus_2 = sfrData(sds_id_2,sds_start, sds_stride, sds_edge, dum2) + Istatus_2
         
      result_real = (( 1.0 - time_wgt) * dum1 ) + (time_wgt * dum2 )
     
      deallocate ( dum1, dum2 )
  
      end function read_nwp_2d 
      
   !  =============================================
   !
   !  ==============================================
      
   function read_nwp_3d (sd_id_1 , sd_id_2, sds_name, time_wgt)  &
      &   result (result_real)
      implicit none
 
         integer, intent(in):: sd_id_1, sd_id_2
      character(len=*), intent(in):: sds_name
	      real, intent(in) :: time_wgt
      integer, dimension(3) :: sds_start, sds_stride, sds_edge
	      real, dimension(:,:,:), allocatable :: result_real
  
      integer :: sds_id_1, sds_id_2
	      real  , dimension(:,:,:), allocatable :: dum1, dum2, dum3
      integer :: istatus_1 = 0
         integer :: istatus_2 = 0
         integer :: istatus= 0
      integer :: gfs_n_lon , gfs_n_lat,  gfs_n_levels
 
      ! - hdf functions
      integer :: sfselect ,  sfn2index, sdreaddata
      integer :: sffattr, sfrnatt
      character (3) ::gfs_array_order_1 , gfs_array_order_2 
	      real  :: MISSING= 9.999E+20
      integer :: i1
  
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF LATITUDES"),gfs_n_lat) + istatus
         istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF LONGITUDES"),gfs_n_lon) + istatus
      istatus = sfrnatt(sd_id_1,sffattr(sd_id_1,"NUMBER OF PRESSURE LEVELS"),gfs_n_levels) + istatus
      istatus = sfrcatt(sd_id_1,sffattr(sd_id_1,"3D ARRAY ORDER"), gfs_array_order_1) + istatus
  
      sds_start = [ 0, 0, 0 ]
      sds_stride = [ 1, 1, 1 ]
  
      if ( gfs_array_order_1 .eq. 'ZXY' ) then
         sds_edge = [ gfs_n_levels, gfs_n_lon , gfs_n_lat ]
         allocate( dum1  (  gfs_n_levels ,gfs_n_lon , gfs_n_lat  ))
         allocate (dum2  (  gfs_n_levels ,gfs_n_lon , gfs_n_lat  ))
         allocate(dum3 (  gfs_n_levels ,gfs_n_lon , gfs_n_lat  ) )
      else 
         sds_edge = [  gfs_n_lon , gfs_n_lat , gfs_n_levels ]
         allocate( dum1  (  gfs_n_lon , gfs_n_lat , gfs_n_levels  ))
         allocate (dum2  (  gfs_n_lon , gfs_n_lat , gfs_n_levels  ))
         allocate(dum3  (  gfs_n_lon , gfs_n_lat , gfs_n_levels   ) )  
      end if
  
      Istatus_1 = 0
         sds_id_1 = sfselect(sd_id_1, sfn2index( sd_id_1, sds_name ) )
      sds_id_2 = sfselect(sd_id_2, sfn2index( sd_id_2, sds_name ) )
  
      !--- read sds's
 
         Istatus_1 = sfrData(sds_id_1,sds_start, sds_stride, sds_edge, dum1) + Istatus_1
         Istatus_2 = sfrData(sds_id_2,sds_start, sds_stride, sds_edge, dum2) + Istatus_2
      
      where ( dum1 .ne. MISSING .and. dum2 .ne. MISSING )
         dum3 = (( 1.0 - time_wgt) * dum1 ) + (time_wgt * dum2 )
      end where 
 
      allocate ( result_real (gfs_n_lon , gfs_n_lat , gfs_n_levels)  )
 
      if ( gfs_array_order_1 .eq. 'ZXY' ) then
         do i1 = 1,   gfs_n_levels
               result_real ( : , : , i1 ) = dum3 (i1, : , :  )  
            end do
      else 
         result_real = dum3
      end if
  
     deallocate ( dum1, dum2, dum3 )
  
      end function read_nwp_3d 

   end subroutine read_gfs_data
   
   !
   !TODO : spatial interpolation 
   !
   subroutine assign_to_sat_grid ( self,  geo)
      use cx_geo_mod
      class ( nwp_main_type ) :: self
      type ( geo_type) :: geo
      integer :: i,j
      allocate ( self % sgrid ( geo % n_x , geo % n_y) ) 
      
      do i = 1 ,  geo % n_x
         do j = 1,  geo % n_y
          
            self % sgrid (i,j) % t_sfc  = self % tmpsfc (geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
            self % sgrid (i,j) % t_tropopause =self % t_trop (geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
            self % sgrid (i,j) % t_air = self % tmpair(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
            self % sgrid (i,j) % p_sfc = self % psfc(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
            self % sgrid (i,j) % rh_sfc =self % rhsfc(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
            self % sgrid (i,j) % p_msl = self % pmsl(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
            self % sgrid (i,j) % tpw =self % tpw(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
            self % sgrid (i,j) % ozone = self % ozone(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % k_index = self % k_index(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % sc_lwp = self % lwp(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )  !TODO this is sc_lwp and not lwp!!
          !  self % sgrid (i,j) % lwp = self % lwp(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % iwp =self % iwp(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % cwp = self % cwp(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % pc =self % pc(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % cloud_fraction =self % cloud_fraction_satellite(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % ncld_layers = self % ncld_layers(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % cld_type =self % cld_type(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % wind_speed = self % wnd_spd_10m(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
          !  self % sgrid (i,j) % wind_dir = self %  wnd_dir_10m(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
           ! self % sgrid (i,j) % lifting_condensation_lev =self %  Lifting_Condensation_Level_Height(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
           ! self % sgrid (i,j) % convection_condensation_lev =self % Convective_Condensation_Level_Height(geo % idx_nwp_x(i,j) , geo % idx_nwp_y(i,j) )
         end do
      end do
       
   end subroutine assign_to_sat_grid
   
   subroutine deallocate_sat_grid ( self)
      class ( nwp_main_type ) :: self
      
      if (allocated ( self % sgrid )) deallocate ( self % sgrid  ) 
   end subroutine deallocate_sat_grid
   
   subroutine deallocate_all ( self )
      class ( nwp_main_type ) :: self
      
      deallocate (self %psfc )
      deallocate (self %lon)
      deallocate (self %lat)

      deallocate (self %pmsl)
      deallocate (self %tmpsfc)
      deallocate (self %tmpair)
      deallocate (self %zsfc)
      deallocate (self %rhsfc)
      deallocate (self %land)
      deallocate (self %weasd)
      deallocate (self %ice)
      deallocate (self %t_trop)
      deallocate (self %p_trop)
      deallocate (self %ozone)
      deallocate (self %tpw)
      
      deallocate ( self%p_std)
      deallocate ( self %t_prof)
      deallocate ( self %z_prof)
      !deallocate ( self %rh_prof)
      deallocate ( self %ozone_prof)
      deallocate ( self %wvmr_prof)
      deallocate (self%wvmrsfc)
      deallocate ( self %sfc_level)
      deallocate ( self % tropo_level)
      
      
   end subroutine deallocate_all
   
   
 
  
end module cx_nwp_mod
