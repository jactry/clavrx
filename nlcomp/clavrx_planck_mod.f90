! $Id:$
!
!
!    All Planck routines 
!   
module clavrx_planck_mod
   private
   public :: planck_tmp2rad
   public :: planck_rad2tmp
      
   integer, parameter:: NPLANCK = 161
   real , parameter :: c1 = 1.191062e-5
   real , parameter :: c2 = 1.4387863
   
   type planck_coef_type
      real :: nu
      real :: a1
      real :: a2
      logical :: is_set      
   end type 
   
   type planck_tbl_type
      real :: b_vec (NPLANCK)
      real :: T_planck ( NPLANCK)
      logical :: is_set
   end type planck_tbl_type
   
   type planck_channel_type
      type ( planck_coef_type ) :: coef
      type ( planck_tbl_type )  :: tbl
      contains
      procedure:: populate_tbl
   end type planck_channel_type
   
   type planck_sensor_type
      logical :: is_set = .false.
      character (len=10):: sensor_set
      type ( planck_channel_type ) :: chn (42) 
      contains
      procedure :: read_coeffs              
   end type
   
   type planck_main_type
      type ( planck_sensor_type ) :: sensor          
   end type
   
   type ( planck_main_type ) :: main
  
   character(len = 20 ) ,save :: sensor_saved 
   
   real, parameter :: T_planck_min = 180.0  
   real, parameter :: delta_T_planck = 1.0
        
contains

   !
   !
   !
   subroutine populate_tbl ( self, sensor )
      implicit none
      class ( planck_channel_type ) :: self
      character ( len = * ) :: sensor
      real :: c1_times_nu__3
      real :: c2_times_nu
            
      integer :: i
      
      if (self % tbl % is_set) return
      
      call main % sensor % read_coeffs ( sensor )
      c1_times_nu__3 =  c1 * self%coef%nu ** 3
      c2_times_nu = c2 * self%coef%nu
      
      do i = 1 , nplanck 
         self % tbl % t_planck(i) = T_planck_min + ( i - 1 ) * delta_T_planck
         self % tbl % B_vec(i) = c1_times_nu__3  / ( exp ( ( c2_times_nu ) / &
              (( self % tbl % T_planck(i) - self%coef%a1 ) / self%coef%a2 ) ) - 1.0) 
      end do
      self % tbl % is_set = .true.
     
   end subroutine populate_tbl 
   
   !
   !
   !
   subroutine read_coeffs ( self , sensor )
      class ( planck_sensor_type ) :: self
      character ( len = * ) :: sensor
      character (len =1024) :: file
      integer :: unit
      character(len = 10 ) :: name
      integer :: chn_idx
      real :: a1 , a2 , nu 
      
      
      file = '/DATA/Ancil_Data/clavrx_ancil_data/avhrr_data/planck_coeff.dat'
      
      
      
      unit = 20
      open ( unit , file = trim(file), action = "read" )
      
      do 
         read(unit, * , end =1 ) name , chn_idx,  a1, a2, nu
        
         if ( trim(name) == trim(sensor)) then   
            self % chn(chn_idx) % coef % a1 = a1
            self % chn(chn_idx) % coef % a2 = a2
            self % chn(chn_idx) % coef % nu = nu
            self % chn(chn_idx) % coef % is_set = .true.
            
         end if   
      end do
      
      
      1 close (unit)
   
   
   end subroutine read_coeffs
   
   !
   !
   !
   function planck_tmp2rad ( tmp , sensor, chn_idx) &
      result ( rad )
      real, intent(in) :: tmp
      character ( len = * ) , intent(in) :: sensor
      integer , intent(in) :: chn_idx
      integer :: l
      real :: b_vec ( NPLANCK) 
      real :: t_vec ( NPLANCK)
      real :: dB_dT_tmp
      
      !
     
      
      
      call main % sensor % chn(chn_idx) % populate_tbl(sensor)
      
      
      b_vec = main % sensor % chn ( chn_idx ) % tbl % b_vec
      t_vec = main % sensor % chn ( chn_idx ) % tbl % t_planck
      l = ( tmp - T_planck_min ) / delta_T_planck
      l = max (1, min ( nplanck - 1 , l ) )
	  
      dB_dT_tmp = (B_vec(l+1)-B_vec(l))/(T_vec(l+1)-T_vec(l))
      rad = b_vec(l) + (tmp - T_vec(l)) * (dB_dT_tmp)
      
   end function planck_tmp2rad
   
   !
   !
   function planck_rad2tmp( rad , sensor , chn_idx ) &
         result ( tmp)
      real, intent(in) :: rad
      character ( len = * ) , intent(in) :: sensor
      integer , intent(in) :: chn_idx
      integer :: l
      real :: b_vec ( NPLANCK) 
      real :: t_vec ( NPLANCK)
      real :: dB_dT_tmp
      
      !
      
      call main % sensor % chn(chn_idx) % populate_tbl(sensor)
      
      b_vec = main % sensor % chn ( chn_idx ) % tbl % b_vec
      t_vec = main % sensor % chn ( chn_idx ) % tbl % t_planck
      call locate ( b_vec, NPLANCK , rad, l )
      l = max (1, min ( nplanck - 1 , l ) )
      
      dB_dT_tmp = (B_vec(l+1)-B_vec(l))/(T_vec(l+1)-T_vec(l))
      tmp = t_vec(l) + ( rad - b_vec ( l)) / dB_dT_tmp 
   
   end function planck_rad2tmp
   
   
   !-------------------------------------------------------------------------
   ! subroutine LOCATE(xx, n, x, j)
   ! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
   !--------------------------------------------------------------------------
   subroutine LOCATE(xx, n, x, j)

   !   Arguments
      integer,                        intent(in)  :: n
      integer,                        intent(out) :: j
      real ,               intent(in)  :: x
      real , dimension(:), intent(in)  :: xx

   !   Local variables
      integer :: i, jl, jm, ju

      jl = 0
      ju = n + 1
      do i = 1, 2*n
         if (ju-jl <= 1) then
            exit
         end if
         jm = (ju + jl) / 2
         if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
      end do
      if (x == xx(1)) then
         j=1
      else if (x == xx(n)) then
         j = n - 1
      else
         j = jl
      end if

   end subroutine LOCATE

end module clavrx_planck_mod

!program test_it

!use clavrx_planck_mod


!real :: p
!p=planck_tmp2rad(233., 'VIIRS',32)

!print*,p

!t= planck_rad2tmp( p, 'VIIRS',32)
!print*,t

!p=planck_tmp2rad(233., 'VIIRS',31)

!print*,p

!t= planck_rad2tmp( p, 'VIIRS',31)
!print*,t

!end
