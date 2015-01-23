! $Header: https://svn.ssec.wisc.edu/repos/aw_clavrx/trunk/imp_mod.f90 17 2014-10-23 06:06:45Z awalther $
!
!  
!   intermediate products
!
module cx_imp_mod
   
   use cx_tools_array_mod,only:array_stats_type,array_statistics,local_radiative_center
   use cx_tools_science_mod
   
   type imp_main_type
      real, allocatable, dimension(:,:) :: bt11_lrc
      real, allocatable, dimension(:,:) :: emis_ch20
      real, allocatable, dimension(:,:) :: refl_sol_ch20
      type (array_stats_type)  :: emis_ch20_3x3
      type (array_stats_type) :: bt_stats_3x3 ( 42)
      type (array_stats_type) :: refl_stats_3x3_ch1
      logical , allocatable :: glint (:,:)
      logical , allocatable :: has_solar_contamination (:,:)
      logical , allocatable :: is_space  (:,:)
      logical , allocatable :: is_valid  (:,:)
   contains
      procedure :: deallocate_all 
      procedure :: populate => populate_imp    
   end type 
   
contains
   !  ------------------------------------------------
   !
   ! -------------------------------------------------
   subroutine deallocate_all ( this)
      class ( imp_main_type) :: this
      if ( allocated ( this % bt11_lrc) ) deallocate ( this % bt11_lrc )
      if ( allocated ( this % emis_ch20) ) deallocate ( this % emis_ch20 )
      if ( allocated ( this % refl_sol_ch20) ) deallocate ( this % refl_sol_ch20 )
      if ( allocated ( this % glint ) ) deallocate ( this % glint)
      if ( allocated ( this % has_solar_contamination ) ) &
               &  deallocate ( this % has_solar_contamination)
      if ( allocated ( this % is_space ) ) deallocate ( this % is_space)
      if ( allocated ( this % is_valid ) ) deallocate ( this % is_valid)
   
   end subroutine deallocate_all
   !  ------------------------------------------------
   !
   ! -------------------------------------------------
   subroutine populate_imp ( this, sat , sfc )
      use cx_sat_mod
      use cx_sfc_mod
      
      type(sat_main_type) , intent(in) :: sat
      type(sfc_main_type) , intent(in) :: sfc
      class (imp_main_type) , intent(in out):: this
      real :: min_bt11 , max_bt11
      integer , parameter :: N_BOX = 3
      character ( len = 10) :: sensor
      real :: sun_earth_distance = 1.
      
      real, parameter :: GLINT_ZEN_THRESH = 40.0
      real, parameter :: FREEZING_TEMP = 273.15
      
      allocate ( this % bt11_lrc (sat % num_pix , sat%num_elem_this_seg) )
     
      call local_radiative_center ( sat % chn(31) % bt  , min_bt11 &
            & , max_bt11 , this % bt11_lrc )
            
      call array_statistics (  sat % chn(31) % bt , n_box , this % bt_stats_3x3(31) )  
      call array_statistics (  sat % chn(1) % ref , n_box , this % refl_stats_3x3_ch1 )  
      ! - emissivity ch20 3.7 um == assumed emissivity computed by ch31 bt
      allocate ( this % emis_ch20 (sat % num_pix , sat%num_elem_this_seg) )
      allocate ( this % refl_sol_ch20 (sat % num_pix , sat%num_elem_this_seg) )
      !sensor ='VIIRS'
      
      
       
      !AW-TODO call compute_emis ( sat % chn(31) % bt, sat % chn(20) % rad , sensor, 20, this % emis_ch20 )
      
      !AW-TODO call compute_solar_refl20 (sat % chn(20) % rad , sat % chn(31) % bt , solar_f0_ch20 &
         !, sat % geo % sol_zen , sun_earth_distance  ,sensor, this % refl_sol_ch20)
      
      call array_statistics ( this % emis_ch20 , n_box, this % emis_ch20_3x3 )
      
      
      
      
      allocate ( this % is_space (sat % num_pix , sat%num_elem_this_seg) )
      this % is_space = .false.
      allocate ( this % is_valid (sat % num_pix , sat%num_elem_this_seg) )
      this % is_valid = .true.
      allocate ( this % has_solar_contamination (sat % num_pix , sat%num_elem_this_seg) )
      this % has_solar_contamination = .false.
      allocate ( this % glint (sat % num_pix , sat%num_elem_this_seg) )
      this % glint = .false.
     
      
      where ( sat % geo % glint_zen < GLINT_ZEN_THRESH &
               & .and. sfc % land_class % data /= ET_land_class % LAND &
               & .and. sfc % snow_class % data == ET_snow_class % NO_SNOW )
         this % glint = .true.
      end where  
     
      ! - exclusions ( freezing pixels, non-uniform, and visible dark
      if ( sat % chan_on (31)) then
         where ( sat % chn(31) % bt < FREEZING_TEMP &
            & .or. this % bt_stats_3x3(31) % std > 1.0 )
         
            this % glint = .false.
         end where
      end if
      
      if ( sat % chan_on (1)) then
         where (  sat % chn(1) % ref < 5.0 )
         
            this % glint = .false.
         end where
      end if 
        
   end subroutine  populate_imp 
   


end module cx_imp_mod
