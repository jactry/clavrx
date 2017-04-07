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
   inp % rfl(2) = 0.28
   inp % rfl(3) = 0.28
   inp % rfl(4) = 0.32
   inp % rfl(5) = 0.33
   inp % rfl(6) = 0.31
   
   inp % sol = 22.
   inp % sat = 13.
   inp % azi = 120.
   
   call inp % info  
   
   call muri_algorithm (inp,out)
   
   print*,out


end program one_pixel_run
