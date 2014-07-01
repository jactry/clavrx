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
      
      cloud_shadow = .false.
      allocate ( distance_km, source = cloud_height) 
      allocate (    Lon_Spacing_Per_m, source = cloud_height)
      
      distance_km = cloud_height * tan (solar_zenith * PI/180.)
      
      Lon_Spacing_Per_m = LAT_SPACING_PER_M * cos ( Lat_pc * PI/180. )
      
      
         do i = 1, 100
            do j = 1,100
               if ( cloud_height (i,j) <= 0. ) cycle
               print*,'===================  ',i,j,'  ===================='   
               print*,distance_km(i,j), solar_zenith (i,j), cloud_height(i,j)
              
               Delta_Lon = sin(Solar_azi(i,j)*Dtor)*distance_km(i,j) * Lon_Spacing_Per_m (i,j)
               Delta_Lat = cos(Solar_azi(i,j)*Dtor)*distance_km (i,j) * Lat_Spacing_Per_m
               print*,delta_lon,delta_lat
               print*,lat(i,j),lat_pc(i,j),lat_pc(i,j) + delta_lat
               print*,lon(i,j),lon_pc(i,j),lon_pc(i,j) + delta_lon
               !print*,lon(i-5:i+5,j-5:j+5)
               !print*,lat(i-5:i+5,j-5:j+5)
            end do
         end do   
      stop
   
   end subroutine cloud_shadow_retr
   
   
   function latlon2pixel ( lat0,lon0,lat , lon , i, j )
      real, intent(in) :: lat0
      real, intent(in) :: lon0
      real, intent(in) :: lat(:,:)
      real, intent(in) :: lon(:,:)
      integer, intent(out) :: i
      integer, intent(out) :: j
      
      
   
   end function latlon2pixel
   

end module cloud_shadow_mod



