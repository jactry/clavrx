!
!
!
module cx_muri_clavrx_bridge_mod

    use muri_interface_mod , only: &
         muri_in_array_type &
       , muri_out_array_type

   
   type cx_muri_structure
      real,allocatable :: aod(:,:)
        contains
        procedure :: allocate =>  cx_muri_structure__allocate
        procedure :: deallocate =>  cx_muri_structure__deallocate
   end type cx_muri_structure
   
   type(cx_muri_structure) :: muri
   
   
contains

   subroutine cx_muri_algorithm 
      implicit none
      type(muri_in_array_type) :: input
      type(muri_out_array_type) :: output
      
      print*,'muri starts'
      
      call  muri_array_loop (input, output )
      muri % aod = output % aod
   end subroutine cx_muri_algorithm
   
   
   
   subroutine cx_muri_structure__allocate(this, dim1, dim2)
      class(cx_muri_structure) :: this
      integer, intent(in) :: dim1,dim2
      allocate ( this % aod ( dim1,dim2))
   end subroutine cx_muri_structure__allocate
   
   subroutine cx_muri_structure__deallocate(this)
      class(cx_muri_structure) :: this
   
   end subroutine cx_muri_structure__deallocate

end module cx_muri_clavrx_bridge_mod
