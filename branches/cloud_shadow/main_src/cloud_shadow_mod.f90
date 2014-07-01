! $Header:$
!
module cloud_shadow_mod


contains

   subroutine cloud_shadow_retr (  &
        cloud_height &
      , solar_azi &
      , solar_zenith &
      , lat &
      , lon &
      , lat_pc &
      , lon_pc &
      , cloud_shadow )
      
      implicit none 
      real, intent(in) :: cloud_height (:,:)
      real, intent(in) :: solar_azi (:,:)
      real, intent(in) :: solar_zenith (:,:)
      real, intent(in) :: lat (:,:)
      real, intent(in) :: lon (:,:)  
      real, intent(in) :: lat_pc (:,:)
      real, intent(in) :: lon_pc (:,:)
      
      logical, intent(out) :: cloud_shadow (:,:)
      
      
      real, parameter :: PI = 3.1415926535897
      real, parameter :: DTOR = PI/180.
      real, allocatable :: distance_km (:,:)
      real, allocatable :: Lon_Spacing_Per_m (:,:)
      real,parameter:: LAT_SPACING_PER_M = 8.9932e-06   ! ( = 1.0/111000.0 m )
      
      integer :: i,j
      real :: delta_lon, delta_lat
      integer :: i_dim, j_dim
      
     
      
      cloud_shadow = .false.
      allocate ( distance_km, source = cloud_height) 
      allocate (    Lon_Spacing_Per_m, source = cloud_height)
      
      distance_km = cloud_height * tan (solar_zenith * PI/180.)
      
      Lon_Spacing_Per_m = LAT_SPACING_PER_M / cos ( Lat_pc * PI/180. )
      
      i_dim = size ( cloud_height, dim=1)
      j_dim = size ( cloud_height, dim=2)
      print*,i_dim,j_dim
         do i = 2 , i_dim - 1
            do j = 2 ,j_dim - 1
               if ( cloud_height (i,j) <= 0. ) cycle
               
               ! are there clear pixels around at all?
               if ( count ( cloud_height (i-1:i+1,j-1:j+1) > 0. ) == 9 ) cycle

               Delta_Lon = sin(Solar_azi(i,j)*Dtor)*distance_km(i,j) * Lon_Spacing_Per_m (i,j)
               Delta_Lat = cos(Solar_azi(i,j)*Dtor)*distance_km (i,j) * Lat_Spacing_Per_m
                     
               call shadow_ind ( lat_pc(i,j) + delta_lat, lon_pc(i,j) + delta_lon , lat, lon , i,j, cloud_shadow) 
               
            end do
         end do   
         
         deallocate ( distance_km)
         deallocate ( lon_spacing_per_m)
      
   
   end subroutine cloud_shadow_retr
   
   
   subroutine shadow_ind ( lat1,lon1, lat , lon , i, j, shad_arr )
      real, intent(in) :: lat1
      real, intent(in) :: lon1
      real, intent(in) :: lat(:,:)
      real, intent(in) :: lon(:,:)
      
      logical, intent(inout) ::shad_arr (:,:)
      integer :: i
      integer :: j
      
      real :: pixel_size_lat (2)
      real :: pixel_size_lon (2)
      integer :: ii,jj
      real :: diff_lat , diff_lon
      real :: delta_lat_ii , delta_lat_jj
      real :: delta_lon_ii , delta_lon_jj
      real :: long_idx,short_idx
      integer :: short_idx_arr
      
   
      
      delta_lat_ii = lat(2,2) - lat(1,2)
      delta_lon_ii = lon(2,2) - lon(1,2)
      delta_lat_jj = lat(2,2) - lat(2,1)
      delta_lon_jj = lon(2,2) - lon(2,1)
      
      diff_lat = lat(i,j) -lat1
      diff_lon = lon(i,j) -lon1
      
      ii = (  diff_lon - ( delta_lon_jj * diff_lat /delta_lat_jj ) ) / &
          ( delta_lon_ii - (delta_lon_jj * delta_lat_ii / delta_lat_jj) )
      
      jj = ( diff_lat - ii * delta_lat_ii) / delta_lat_jj
      
      
      long_idx   = maxval (abs([ii,jj]))
      short_idx  = minval (ABS([ii,jj]))
      
      do k = 1 , long_idx      
         
         
         if (abs(ii) == long_idx ) then
            short_idx_arr = CEILING ( short_idx * sign(k,ii) / long_idx  ) 
            shad_arr(i+ sign(k,ii), j+ short_idx_arr) = .true.
            shad_arr(i+ sign(k,ii), j+ short_idx_arr-1) = .true.
         else 
            short_idx_arr = CEILING ( short_idx * sign(k,jj) / long_idx  ) 
            shad_arr(i+short_idx_arr,j + sign ( k,jj) )= .true.
            shad_arr(i+short_idx_arr-1,j + sign(k,jj) )= .true.
         endif      
      end do
      
      
    
      
   
   end subroutine  shadow_ind
   
   
   
   

end module cloud_shadow_mod



