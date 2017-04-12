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
   
   integer :: i,ii
   
   
   
   do i=20,20
   
   
   
   inp % rfl(1) = 0.3 
   inp % rfl(2) = 0.28
   inp % rfl(3) = 0.28
   inp % rfl(4) = i/100.   !0.32
   inp % rfl(5) = 0.33
   inp % rfl(6) = 0.31
   
   inp % sol = 22.
   inp % sat = 13.
   inp % azi = 120.
   
   inp % ws = 2.4
   
   !call inp % info  
   
   call muri_algorithm (inp,out)
   print*
   
   print*,'channel 4 reflectance: ',inp % rfl(4) 
   print*,'cm mode, fm mode :',out % cm_mode, out % fm_mode, out % fmf
   print*,'AOT reference: ',out % aot
   do ii=1,6 
   
   print*,'AOT channel',ii,out % aot_channel(ii)
   
   end do
   end do
   
   
  


end program one_pixel_run
