! $Header:$
module cloud_type_algo_module
   
   type  et_cloudtype_class_type
      integer :: CLEAR = 0
      integer :: PROB_CLEAR = 1
      integer :: FOG = 2
      integer :: WATER = 3
      integer :: SUPERCOOLED = 4
      integer :: MIXED = 5
      integer :: OPAQUE_ICE = 6
      integer :: TICE = 6
      integer :: CIRRUS = 7
      integer :: OVERLAP = 8
      integer :: OVERSHOOTING = 9
      integer :: UNKNOWN = 10
      integer :: DUST = 11
      integer :: SMOKE = 12
      integer :: FIRE = 13
   end type 
   
   type cloud_type_sat_type
      logical , dimension(42) :: chan_on
      real :: bt_ch20
      real :: bt_ch27
      real :: bt_ch29
      real :: bt_ch31
      real :: bt_ch32
      real :: ref_ch1
      real :: ref_ch2
      real :: ref_ch6
      real :: ref_ch7
      real :: ref_ch8
      real :: ref_ch26
      real :: emis_ch20_3x3_mean
      real :: ref_ch1_3x3_std
      real :: ref_ch1_3x3_min
      
   end type cloud_type_sat_type

   
   type cloud_type_input_type
      type ( cloud_type_sat_type) :: sat
   
   end type cloud_type_input_type

contains
   subroutine cloud_type_algo
      print*,'start cloud type algo'
   
   
   end subroutine cloud_type_algo


end module cloud_type_algo_module
