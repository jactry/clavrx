! $Id$
!
!

subroutine muri_array_loop (input, output )
   
   use muri_interface_mod,only: &
      muri_in_array_type &
      , muri_out_array_type
   
   use muri_retrieval_mod
   
   use muri_definitions_mod, only: &
      muri_input_type  &
      , muri_output_type
   
   implicit none
   
   type ( muri_in_array_type) , intent(in):: input
   type ( muri_out_array_type) :: output
   type ( muri_input_type) :: inp_pixel
   type ( muri_output_type) :: out_pixel
   integer :: i,j,k
   
   print*, '========>    .. >',input % dim
  
   do i = 1, input % dim(1)
  
      do j = 1, input % dim(2)
         if ( .not. input % do_it(i,j)) cycle
        
         inp_pixel % sat = input % sat(i,j)
         inp_pixel % sol = input % sol(i,j)
         inp_pixel % azi = input % azi(i,j)
         inp_pixel % ws = 2.0
         do k=1,6
            inp_pixel % rfl(k) = input % ref(k,i,j)
         end do
        ! call inp_pixel % info
         call muri_algorithm( inp_pixel, out_pixel )
      !   print*,'done ..'
       !  do k=1, 6
           ! print*,out_pixel % aot
            !output % aot_channel(k,i,j) = out_pixel % aot_channel(k)
       !  end do 
       !  print*,'w'
         output % aot(i,j) = out_pixel % aot
        ! output % cm_mode(i,j) = out_pixel % cm_mode
         !output % fm_mode(i,j) = out_pixel % fm_mode
        ! print*,'w3'
        ! output % fmf(i,j) = out_pixel % fmf 
      end do
   end do   
   
  
   
  

   




end subroutine muri_array_loop
