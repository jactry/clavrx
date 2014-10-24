! $Id: pixel_routines.f90 538 2014-09-15 20:30:48Z dbotambekov $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: acha_shadow_module.f90 (src)
!       ACHA_SHADOW (program)
!
! PURPOSE: this module computes the geometrical shadow routine
!
! DESCRIPTION: 
!
! AUTHORS:  Andi Walther, UW/SSEC/CIMSS
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! Public routines used in this MODULE:
!
!--------------------------------------------------------------------------------------
module ACHA_SHADOW

 implicit none

 public:: CLOUD_SHADOW_RETR

 private:: SHADOW_IND          

 contains

!==============================================================================
! Cloud Shadow Routine
!
! -1 = shadow estimation not possible
!  0 = no shadow
!  1 = shadow determined from geometrical analysis
!  2 = shadow determined from spectral analysis
!  3 = shadow determined from both geometrical and spectral analysis
!
!==============================================================================
subroutine CLOUD_SHADOW_RETR (  &
           Cloud_Height &
         , Solar_Azi &
         , Solar_Zenith &
         , Lat &
         , Lon &
         , Lat_Pc &
         , Lon_Pc &
         , Cloud_Shadow )
      
      implicit none 
      real, intent(in) :: Cloud_Height (:,:)
      real, intent(in) :: Solar_Azi (:,:)
      real, intent(in) :: Solar_Zenith (:,:)
      real, intent(in) :: Lat (:,:)
      real, intent(in) :: Lon (:,:)  
      real, intent(in) :: Lat_Pc (:,:)
      real, intent(in) :: Lon_Pc (:,:)
      
      integer(1), intent(out) :: Cloud_Shadow (:,:)
      
      real, parameter :: PI = 3.1415926535897
      real, parameter :: DTOR = PI/180.
      real, allocatable :: Distance_km (:,:)
      real, allocatable :: Lon_Spacing_Per_m (:,:)
      real,parameter:: LAT_SPACING_PER_M = 8.9932e-06   ! ( = 1.0/111000.0 m )
      
      integer :: i,j
      real :: delta_lon, delta_Lat
      integer :: i_dim, j_dim
      
      !--- initialize output
      cloud_shadow = 0

      !--- allocate local memory
      i_dim = size ( Cloud_Height, dim=1)
      j_dim = size ( Cloud_Height, dim=2)
      allocate ( Distance_km (i_dim, j_dim))
      allocate ( Lon_Spacing_Per_m (i_dim, j_dim))

      Distance_km = Cloud_Height * tan (Solar_Zenith * DTOR )
      
      Lon_Spacing_Per_m = LAT_SPACING_PER_M / cos ( Lat_Pc * DTOR )
      
         do i = 2 , i_dim - 1
            do j = 2 ,j_dim - 1
               if ( Solar_Zenith (i,j) > 85.0 ) cycle
               if ( Cloud_Height (i,j) <= 0.0 ) cycle
               
               ! are there clear pixels around at all?
               if ( count ( Cloud_Height (i-1:i+1,j-1:j+1) > 0.0 ) == 9 ) cycle

               Delta_Lon = -1.0 * sin(Solar_Azi(i,j) * DTOR ) * Distance_km(i,j) * Lon_Spacing_Per_m(i,j)
               Delta_Lat = -1.0 * cos(Solar_Azi(i,j) * DTOR ) * Distance_km(i,j) * Lat_Spacing_Per_m
                
               call SHADOW_IND ( Lat_Pc(i,j) + Delta_Lat, Lon_Pc(i,j) + Delta_Lon, Lat, Lon, i, j, Cloud_Shadow) 
              
            end do
         end do   
         
         !--- destroy local memory
         deallocate ( Distance_km)
         deallocate ( Lon_Spacing_Per_m)

   end subroutine CLOUD_SHADOW_RETR
   !================================================================================================================ 
   !  
   !  input:
   !     Lat1,lon1 : the lon/Lat values of maximal shadow
   !
   !     Lat , Lon : lon/Lat coordinate array
   !    
   !      i , j : the cloud pixel for Lat/Lon array
   !
   !  output : 
   !      shad_arr: one byte integer array with same size as Lat, lon (0 = false, 1 = true)
   !
   !================================================================================================================ 
   subroutine SHADOW_IND( Lat1,lon1, Lat , Lon , i, j, shad_arr )
      real, intent(in) :: Lat1
      real, intent(in) :: lon1
      real, intent(in) :: Lat(:,:)
      real, intent(in) :: Lon(:,:)
      
      integer(1) , intent(inout) ::shad_arr (:,:)
      integer , intent(in) :: i
      integer , intent(in) :: j
      
      real :: pixel_size_Lat (2)
      real :: pixel_size_Lon (2)
      integer :: ii,jj
      real :: diff_Lat , diff_lon
      real :: delta_Lat_ii , delta_Lat_jj
      real :: delta_lon_ii , delta_lon_jj
      integer :: long_idx,short_idx
      integer :: short_idx_arr
      
      integer :: dim_1 , dim_2
      integer :: idx_1 , idx_2
      integer :: k
      
      dim_1 = size ( shad_arr, 1 )
      dim_2 = size ( shad_arr, 2 )
      
      delta_Lat_ii = Lat(i,j) - Lat(i-1,j)
      delta_lon_ii = Lon(i,j) - Lon(i-1,j)
      delta_Lat_jj = Lat(i,j) - Lat(i,j-1)
      delta_lon_jj = Lon(i,j) - Lon(i,j-1)
     
      diff_Lat = Lat1 - Lat(i,j)
      diff_Lon = lon1 - Lon(i,j) 
      
      !     use of these equations:
      ! diff_Lon = ii * delta_lon_ii + jj * delta_lon_jj
      ! diff_Lat = ii * delta_Lat_ii + jj * delta_Lat_jj
      !   solve for ii => 
      
      ii = ( diff_Lat * delta_lon_jj - diff_Lon * delta_Lat_jj) / &
         (delta_Lat_ii * delta_lon_jj - delta_lon_ii * delta_Lat_jj)
      
      jj = ( diff_Lat - ii * delta_Lat_ii) / delta_Lat_jj
     
      ! - find longer dim
      !
      long_idx   = maxval (ABS([ii,jj]))
      short_idx  = minval (ABS([ii,jj]))
      
      ! - this fills all pixel from Lat1 to Lat(i,j)/Lon(i,j)  
      do k = 1 , long_idx        
         if (abs(ii) == long_idx ) then
            idx_1 = i + sign(k,ii)
            if (idx_1 < 1 .or. idx_1 > dim_1 ) cycle
             
            short_idx_arr = CEILING ( short_idx * sign(k,jj) / (1.* long_idx ) ) 
       
            idx_2 = j + short_idx_arr
            if ( idx_2 > 0 .and. idx_2 <= dim_2 ) shad_arr( idx_1 , idx_2 )   = 1
            idx_2 = j + short_idx_arr - 1
            if ( idx_2 > 0 .and. idx_2 <= dim_2 ) shad_arr( idx_1 , j + short_idx_arr-1) = 1
         else 
            idx_2 = j + sign(k,jj)
          
            if (idx_2 < 1 .or. idx_2 > dim_2 ) cycle
            short_idx_arr = CEILING ( short_idx * sign(k,ii) / ( 1.* long_idx)  ) 
            
            idx_1 = i + short_idx_arr
            
            if ( idx_1 > 0 .and. idx_1 <= dim_1 ) shad_arr(i+short_idx_arr,idx_2 ) = 1
            idx_1 = i + short_idx_arr - 1
            if ( idx_1 > 0 .and. idx_1 <= dim_1 ) shad_arr(i+short_idx_arr-1, idx_2) = 1
         endif      
      end do
   
   end subroutine SHADOW_IND

!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module ACHA_SHADOW
