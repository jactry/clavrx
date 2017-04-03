! $Id$

module muri_retrieval_mod
   use muri_definitions_mod, only: &
      muri_input_type &
      , muri_output_type
   
   use muri_forward_mod, only: &
      muri_fwd_type &
      , muri_forward
      
   implicit none
   private
   public :: muri_algorithm 
contains
   subroutine muri_algorithm (inp, out)
      type( muri_input_type), intent(in) :: inp
      type( muri_output_type), intent(out) :: out
      
      type ( muri_fwd_type ) :: fwd
      
      
      print*,'start muri algorithm'
      call inp%info
      
      
      ! - first prepare input
      ! - - done
      
      
      ! - now prepare the table
      
      
      
      call muri_forward ( 4.5, fwd)
      
      
      ! - do something with forward model
      
      
      
      
       !- once find a good match between fwd and measurment give output
      out% aot = 3.4
      
      
      
   
   
   end subroutine  muri_algorithm  



end module muri_retrieval_mod
