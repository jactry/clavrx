! $Id$
!
!
module muri_interface_mod

   type muri_in_array_type
      integer :: dim(2)
      real,allocatable :: sol(:,:)
      real,allocatable :: sat(:,:)
      real,allocatable :: azi(:,:)
      real,allocatable :: ref(:,:,:)
      real,allocatable :: windspeed(:,:)
      
      contains
      procedure :: allocate=>muri_in_array_type__allocate
      procedure :: deallocate=>muri_in_array_type__deallocate
   
   end type muri_in_array_type
   
   type muri_out_array_type
      real :: aod
   
   end type muri_out_array_type

contains

   subroutine muri_in_array_type__allocate(this,dim1,dim2)
      class(muri_in_array_type) :: this
      integer :: dim1,dim2
      this % dim(1) = dim1
      this % dim(2) = dim2
      allocate ( this % sol( dim1,dim2))
      allocate ( this % sat( dim1,dim2))
      allocate ( this % azi( dim1,dim2))
      allocate ( this % windspeed( dim1,dim2))
      
     
      allocate ( this % ref( 6,dim1,dim2))
      
   
   end subroutine  muri_in_array_type__allocate  
   
   
   subroutine muri_in_array_type__deallocate(this)
      class(muri_in_array_type) :: this
      this % dim = [0,0]
      deallocate ( this % sol)
      deallocate ( this % sat)
      deallocate ( this % azi)
      deallocate ( this % windspeed)
      
     
      deallocate ( this % ref)
      
   
   end subroutine  muri_in_array_type__deallocate  

end module muri_interface_mod
