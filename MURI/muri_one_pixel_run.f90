! $Id$
program one_pixel_run

   use muri_retrieval_mod
   
   use muri_definitions_mod, only: &
      muri_input_type  &
      , muri_output_type
   
   ! ---------------------------------------   
      
   implicit none   
   type ( muri_input_type) :: inp
   type ( muri_output_type) :: out
   inp % rfl(1) = 0.3   
   
   call muri_algorithm (inp,out)
   
   print*,out


end program one_pixel_run
