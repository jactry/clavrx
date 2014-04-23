! $Header:$


module cloud_type_bridge_module

   use CONSTANTS
   use PICEL_COMMON
   use CLOUD_TYPE_ALGO_MODULE 
   
   implicit none
   
   public :: CLOUD_TYPE_BRIDGE  

contains
   subroutine cloud_type_bridge
      implicit none
      print*,'start type bridge' 
      
      call cloud_type_algo
      
      stop
   
   end subroutine cloud_type_bridge

end module cloud_type_bridge_module
