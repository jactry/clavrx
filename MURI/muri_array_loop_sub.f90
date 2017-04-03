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
   print*,'starrrrrt...'
   print*,shape(input)
   
   
  do i=1, 6
   print*,input % ref (i,2300,50)
  end do
  
  
  do i = 1, input % dim(1)
      do j = 1, input % dim(2)
         inp_pixel % sat = input % sat(i,j)
         inp_pixel % sol = input % sol(i,j)
         inp_pixel % azi = input % azi(i,j)
         do k=1,6
            inp_pixel % rfl(k) = input % ref(k,i,j)
         end do
         call muri_algorithm( inp_pixel, out_pixel )
      
      end do
   end do   
   
   print*,'end....'






end subroutine muri_array_loop
