!
! $Id:$
!
!   Usage:
!
!  
!
!
module cx_tools_science_planck_mod
   implicit none
      
   private
  
   public :: planck_bt2rad
   public :: planck_rad2bt
   
   ! -  communication
   character (len = 9 ) , save :: sensor_saved 
   
   !- constants
   integer , parameter :: N_PLANCK = 161
   real :: nu_20 
   
   type :: chn_type
      logical :: is_set  
      real :: nu
      real :: a1
      real :: a2
      real, dimension (N_PLANCK ) :: rad_table
      real, dimension (N_PLANCK ) :: tmp_table
   end type chn_type
   
   type ( chn_type ), dimension (20:39) , save :: coef
   
   interface planck_bt2rad
      module procedure &
          planck_rad_2d &
         , planck_rad_skalar
      
   
   end interface planck_bt2rad
   
     interface planck_rad2bt
      module procedure &
          planck_temp_2d &
         , planck_temp_skalar
      
   
   end interface planck_rad2bt
   
contains
   subroutine read_coef (sensor)
      character(len=*), intent(in) :: sensor
      
      character ( len =200) :: a
      character ( len =200) :: b
      integer :: chn ,  i
      real, parameter :: t_planck_min  = 180.
      real, parameter :: delta_t_planck = 1.0
      real, parameter :: c1 = 1.191062e-05
      real, parameter :: c2 = 1.4387863
      integer :: j_chn
      
      coef(:)%is_set = .false.
!TODO -cool make a generic tool from this!!!      
      open(unit=10,file='planck_coef.dat')
      do while ( 1 == 1)
         read(unit=10,end = 200 , fmt=*) a
       
         if ('S-'//trim(sensor) == trim(a)) then
            do while ( b(:1) /='S') 
               read(unit=10, end = 200 , fmt=*) b
               if ( b(:1) /= 'S') then
                 
                  read (b(4:5),*) chn 
                     
                  read(unit=10, fmt = * ) coef (chn) % a1 ,  coef (chn) % a2 ,  coef (chn) % nu
                  coef(chn)% is_set = .true.
               end if
            end do   
         end if
      
      end do
      
     
      200 continue
      do j_chn = 20, 39 
         if ( .not. coef(j_chn)%is_set ) cycle
         do i = 1, N_PLANCK
            coef (j_chn) % tmp_table (i) = t_planck_min + ( i -1 ) * delta_t_planck
            coef ( j_chn) % rad_table (i) = c1 * (coef (j_chn) % nu ** 3) / ( exp ( ( c2 * coef(j_chn) % nu ) / &
               (( coef (j_chn) % tmp_table(i) - coef(j_chn) % a1) / coef(j_chn) % a2 )) -1.0 )
              
         end do
      end do
      
     
      
      sensor_saved = trim(sensor)
     
   
   end subroutine read_coef
   !
   !
   !
   subroutine planck_rad_2d ( bt , sensor, chn , rad ) 
      real ,dimension(:,:) , intent(in) :: bt
      character ( len = *) , intent (in) :: sensor
      integer , intent(in) :: chn
      
      real, dimension (:,:) , intent(out) :: rad
      real, dimension(n_planck) :: T_planck
      real, parameter :: T_planck_min = 180.0
      real, parameter :: delta_T_planck = 1.0
      integer :: i 
      integer :: l
      
      integer, dimension(2) :: dim 
      integer :: k , m
      real :: dB_dT_tmp
    
      if ( trim(sensor) /= trim(sensor_saved) ) then
         call read_coef ( sensor )
      end if
          
   	  
     
      dim = shape ( bt)
      do k =1 , dim(1)
         do m = 1, dim(2)
            l = ( bt (k,m) - T_planck_min ) / delta_T_planck
            l = max (1, min ( n_planck - 1 , l ) )
            dB_dT_tmp = ( coef(chn) % rad_table (l +1) - coef(chn) % rad_table ( l ) ) &
                      /  delta_t_planck
	  
            rad ( k , m)  = coef(chn) % rad_table( l ) + ( bt( k , m ) - coef( chn ) % tmp_table (l) ) * (dB_dT_tmp)
         end do         
      end do
      

   end subroutine planck_rad_2d
   !
   !
   !
   subroutine planck_rad_skalar ( bt , sensor, chn , rad ) 
      real  , intent(in) :: bt
      character ( len = *) , intent (in) :: sensor
      integer , intent(in) :: chn
      
      real,  intent(out) :: rad
      real, dimension(n_planck) :: T_planck
      real, parameter :: T_planck_min = 180.0
      real, parameter :: delta_T_planck = 1.0
      integer :: i 
      integer :: l
      
      
      
      real :: dB_dT_tmp
    
      if ( trim(sensor) /= trim(sensor_saved) ) then
         call read_coef ( sensor )
      end if
          
   	  
     
     
      
            l = ( bt - T_planck_min ) / delta_T_planck
            l = max (1, min ( n_planck - 1 , l ) )
            dB_dT_tmp = ( coef(chn) % rad_table (l +1) - coef(chn) % rad_table ( l ) ) &
                      /  delta_t_planck
	  
            rad  = coef(chn) % rad_table( l ) + ( bt- coef( chn ) % tmp_table (l) ) * (dB_dT_tmp)
         
      

   end subroutine planck_rad_skalar
   
   
   
   !
   !
   !
   subroutine planck_temp_2d ( rad , sensor, chn , bt ) 
      real ,dimension(:,:) , intent(in) :: rad
      character ( len = *) , intent (in) :: sensor
      integer , intent(in) :: chn
      
      real, dimension (:,:) , intent(out) :: bt
      real, dimension(n_planck) :: T_planck
      real, parameter :: T_planck_min = 180.0
      real, parameter :: delta_T_planck = 1.0
      integer :: i 
      integer :: l
   
      integer, dimension(2) :: dim 
      integer :: k , m
      real :: dB_dT_tmp
      if ( sensor /= sensor_saved ) then
         call read_coef ( sensor )
      end if
   	  
      
      dim = shape ( rad )
      do k =1 , dim(1)
         do m = 1, dim(2)
           
            
            l = minloc(abs( coef(chn) % rad_table - rad(k,m) ),1)
            l = max (1, min ( n_planck - 1 , l ) )
           
            dB_dT_tmp = ( coef(chn) % rad_table (l+1) - coef(chn) % rad_table(l) )  /  delta_t_planck	  
            bt ( k , m)  = coef(chn)%tmp_table(l) + (rad(k,m) - coef(chn)%rad_table(l)) / (dB_dT_tmp)
         end do         
      end do
      

   end subroutine planck_temp_2d
   
    subroutine planck_temp_skalar ( rad , sensor, chn , bt ) 
      real , intent(in) :: rad
      character ( len = *) , intent (in) :: sensor
      integer , intent(in) :: chn
      
      real, intent(out) :: bt
      real, dimension(n_planck) :: T_planck
      real, parameter :: T_planck_min = 180.0
      real, parameter :: delta_T_planck = 1.0
      
      integer :: l
   
     
      real :: dB_dT_tmp
      if ( sensor /= sensor_saved ) then
         call read_coef ( sensor )
      end if
   	  

      l = minloc(abs( coef(chn) % rad_table - rad ),1)
      l = max (1, min ( n_planck - 1 , l ) )
           
      dB_dT_tmp = ( coef(chn) % rad_table (l+1) - coef(chn) % rad_table(l) )  /  delta_t_planck	  
      bt   = coef(chn)%tmp_table(l) + (rad - coef(chn)%rad_table(l)) / (dB_dT_tmp)
 
   end subroutine planck_temp_skalar
   
  
   
end module cx_tools_science_planck_mod


!program test_it
!   use planck_mod
   !-- user requirements
   ! 1. compute bt from rad
   ! 2. compute rad from bt
   ! 3. input should be channel-nr and sensor name
!   real , dimension(2,2) :: mmm
!   real, dimension(2,2) :: nnn
!   real :: m , n
!   nnn(1,:) =[290.1,269.3042]
!   nnn(2,:) = [322.2,311.99]
!   print*,nnn
!   call  planck_bt2rad( nnn , 'VIIRS' , 31 , mmm )
!   print*,mmm
!call planck_rad2bt ( mmm , 'VIIRS' , 31 , nnn )

!n=298.4
!call  planck_bt2rad( n , 'VIIRS' , 31 , m )
! print*,m
!call planck_rad2bt ( m , 'VIIRS' , 31 , n )
!print*,n
!print*,n
!end program test_it
