! $Header:$


module cloud_type_bridge_module

   use CONSTANTS
   use PIXEL_COMMON
   use CLOUD_TYPE_ALGO_MODULE 
   
   implicit none
   
   public :: CLOUD_TYPE_BRIDGE  

contains
   subroutine cloud_type_bridge
      implicit none
      
      type ( cloud_type_input_type) :: type_inp
      integer :: i , j
      
      print*,'start type bridge' 
      
      ! -----------    loop over pixels -----   
      line_loop: do i = 1, num_pix
         elem_loop: do  j = 1,num_scans_read
            call cloud_type_algo
         end do elem_loop
      end do   line_loop
      
      
      
      stop
   
   end subroutine cloud_type_bridge

end module cloud_type_bridge_module
