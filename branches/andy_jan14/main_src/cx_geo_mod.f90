module cx_geo_mod

   type geo_type
      integer , dimension(:,:), allocatable :: idx_nwp_x
      integer , dimension(:,:), allocatable :: idx_nwp_y
      integer :: n_x
      integer :: n_y
      
      
      contains
      procedure :: deallocate_all => deallocate_geo
   
   end type geo_type
   
 

contains

   ! --
   ! ---
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   !   lat_master lon_master: 
   !      coarser grid ( NWP)
   !    lat_slave, lon_slaver:
   !       higher resoluted grid ( SAT)
   !       geo % idx_nwp_x , geo % idx_nwp_y 
   !       master index in slave grid.
   !      
   !
   !
   subroutine lon_lat_index &
         ( lat_master , lon_master &
         , lat_slave  , lon_slave &
         , geo )
      
      implicit none
      
      real, dimension(:,:) , intent(in) :: lat_master , lon_master
      real, dimension(:,:) , intent(in) :: lat_slave  , lon_slave  
      type ( geo_type ) , intent(out) :: geo
   
      integer , dimension(2) :: dim_ma 
      integer , dimension(2) :: dim_sl 
      integer :: dim_ma_1 , dim_ma_2
      integer :: dim_sl_1 , dim_sl_2
   
      real, dimension(:,:), allocatable :: lat_loc
      real, dimension(:,:), allocatable :: lon_loc
       integer :: i,j
      real :: dLat_ma, dLon_ma
  
      dim_ma = shape ( lat_master )
      dim_ma_1 = dim_ma (1)
      dim_ma_2 = dim_ma (2)
      dim_sl = shape ( lat_slave )
      dim_sl_1 = dim_sl (1)
      dim_sl_2 = dim_sl (2)
      
      geo % n_x = dim_sl_1
      geo % n_y = dim_sl_2
		geo % idx_nwp_x = 1
		geo % idx_nwp_y = 1
      allocate ( geo % idx_nwp_x (dim_sl_1 , dim_sl_2 ) ,  geo % idx_nwp_y (dim_sl(1) , dim_sl(2) ))
      allocate ( lat_loc (dim_sl_1 , dim_sl_2 ) ,  lon_loc (dim_sl(1) , dim_sl(2) ))
      
      lat_loc = lat_slave
      lon_loc = lon_slave
      
      where ( lon_loc < 0 )
         lon_loc = lon_loc + 360
      end where
      
      dLat_ma = lat_master(10,11)  - lat_master(10,10)
      dLon_ma = lon_master(11,10) - lon_master(10,10)
		
		
      geo % idx_nwp_x = nint((( lon_loc - lon_master(1,1) ) / dlon_ma ) + 1)
      geo % idx_nwp_y = nint((( lat_loc - lat_master(1,1) ) / dlat_ma ) + 1)
        
      where ( geo % idx_nwp_x > dim_ma_1 )
         geo % idx_nwp_x = 1
      end where       
      
      deallocate ( lat_loc , lon_loc )
      
      ! - TODO This is very simple
      !  - I may have to adjust this ...
      
   end subroutine lon_lat_index
   

   !
   !
   !
    subroutine deallocate_geo (  self )
      class ( geo_type )  :: self
      
      if ( allocated (self % idx_nwp_x)) deallocate ( self % idx_nwp_x)
      if ( allocated (self % idx_nwp_y))  deallocate ( self % idx_nwp_y)
    
    
    
    end subroutine deallocate_geo
   
end module cx_geo_mod
