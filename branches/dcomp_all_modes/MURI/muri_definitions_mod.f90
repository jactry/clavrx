! $Id$

module muri_definitions_mod

    type muri_input_type
      real :: rfl(6) 
      real :: sol
      real :: sat
      real :: azi  
      
      contains
      
      procedure :: info => muri_input_info
      
   end type muri_input_type
   
   type muri_output_type
      real :: aot
   
   end type  muri_output_type
   
   
   
contains

   subroutine muri_input_info (this) 
      class(muri_input_type) :: this
      integer :: i
      print*,'=====   MURI INPUT TYPE'
      print*,'REFELCTANCES:'
      do i=1,6 
         print*,'channel ',i,this% rfl(i)
      end do
      print*,'solar zenith: ',this%sol
      print*,'satellite zenith: ',this%sat
      print*,'relative azimuth difference: ',this%azi
      print*
   
   end subroutine  muri_input_info 
   

end module muri_definitions_mod
