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
      logical, allocatable :: do_it(:,:)
      
      contains
      procedure :: allocate=>muri_in_array_type__allocate
      procedure :: deallocate=>muri_in_array_type__deallocate
   
   end type muri_in_array_type
   
   type muri_out_array_type
      integer :: dim(2)
      real, allocatable :: aot(:,:)
      real, allocatable :: aot_channel(:,:,:)
      real, allocatable :: fmf(:,:)
      integer, allocatable :: fm_mode(:,:)
      integer, allocatable :: cm_mode(:,:)
      
      
      contains
      procedure :: allocate=>muri_out_array_type__allocate
      procedure :: deallocate =>muri_out_array_type__deallocate
   
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
      allocate ( this % do_it(dim1,dim2))
      
     
      allocate ( this % ref( 6,dim1,dim2))
      
   
   end subroutine  muri_in_array_type__allocate  
   
   
   subroutine muri_in_array_type__deallocate(this)
      class(muri_in_array_type) :: this
      this % dim = [0,0]
      deallocate ( this % sol)
      deallocate ( this % sat)
      deallocate ( this % azi)
      deallocate ( this % windspeed)
      deallocate ( this % do_it)
     
      deallocate ( this % ref)
      
   
   end subroutine  muri_in_array_type__deallocate 
   
   
   
   !
   !
   !
   subroutine muri_out_array_type__allocate(this,dim1,dim2)
      class(muri_out_array_type) :: this
      integer, intent(in) :: dim1,dim2
      this % dim(1) = dim1
      this % dim(2) = dim2
      allocate ( this % aot( dim1,dim2))
      allocate ( this % aot_channel(6, dim1,dim2))
      allocate ( this % fmf( dim1,dim2))
      allocate ( this % fm_mode( dim1,dim2))
      allocate ( this % cm_mode( dim1,dim2))
      this % aot = -999.
   
   end  subroutine muri_out_array_type__allocate
   
   
   !
   !
   !
      subroutine muri_out_array_type__deallocate(this)
      class(muri_out_array_type) :: this
      this % dim = [0,0]
      deallocate ( this % aot)
      deallocate ( this % aot_channel)
      deallocate ( this % fmf)
      deallocate ( this % fm_mode)
      deallocate ( this % cm_mode)
     
      
      
   
   end subroutine  muri_out_array_type__deallocate 

end module muri_interface_mod
